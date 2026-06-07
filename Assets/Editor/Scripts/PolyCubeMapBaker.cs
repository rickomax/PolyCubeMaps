// PolyCubeMapBaker.cs
// Editor utility for generating PolyCube-Map data (Tarini et al. 2004) that is
// compatible with the runtime shader at Assets/Resources/PolyCubeMap.shader and
// the MonoBehaviour applier at Assets/Scripts/PolyCubeMap.cs.
//
// The shader expects two things packed into one RGBA32 texture:
//   1) A small "LUT" region whose pixels encode, per polycube vertex cell,
//      (patchX_in_squarelets, patchY_in_squarelets, encodedCaseAndRotation, 255).
//   2) A larger region of square "patches" that contain the actual color data,
//      sliced into S x S squarelets.
//
// We implement (per paper §3, §5, §6.1) only the cases 3, 4a, 4b, 5 -- cases
// 6a / 6b are intentionally skipped (the shader's fragment program does not
// implement their projection branches, just like the reference impl).

using System;
using System.Collections.Generic;
using System.IO;
using UnityEditor;
using UnityEngine;

public static class PolyCubeMapBaker
{
    // ----------------------------------------------------------------------
    // 1) PolyCube data type
    // ----------------------------------------------------------------------

    // A polycube is just a set of axis-aligned unit cubes identified by their
    // integer minimum corner. The polycube *vertices* live on the integer
    // grid (the corners of those cubes); the LUT (per paper §3) is built over
    // the dual partition, i.e. one LUT entry per polycube vertex.
    public class PolyCube
    {
        public List<Vector3Int> Cubes = new List<Vector3Int>();

        public PolyCube() { }
        public PolyCube(IEnumerable<Vector3Int> cubes) { Cubes.AddRange(cubes); }

        // All vertices that bound at least one cube (i.e. potentially non-empty
        // dual cells). For a single unit cube at the origin, this yields 8.
        public HashSet<Vector3Int> EnumerateVertices()
        {
            var set = new HashSet<Vector3Int>();
            foreach (var c in Cubes)
            {
                for (var dz = 0; dz <= 1; dz++)
                for (var dy = 0; dy <= 1; dy++)
                for (var dx = 0; dx <= 1; dx++)
                    set.Add(new Vector3Int(c.x + dx, c.y + dy, c.z + dz));
            }
            return set;
        }

        public Bounds VertexBounds()
        {
            if (Cubes.Count == 0) return new Bounds();
            var minC = Cubes[0];
            var maxC = Cubes[0] + Vector3Int.one;
            foreach (var c in Cubes)
            {
                minC = Vector3Int.Min(minC, c);
                maxC = Vector3Int.Max(maxC, c + Vector3Int.one);
            }
            var b = new Bounds();
            b.SetMinMax((Vector3)minC, (Vector3)maxC);
            return b;
        }

        public bool ContainsCube(Vector3Int c)
        {
            // tiny membership test; polycubes here are small (<32 cubes).
            for (var i = 0; i < Cubes.Count; i++) if (Cubes[i] == c) return true;
            return false;
        }
    }

    // ----------------------------------------------------------------------
    // 2) Cell-config + rotation lookup
    // ----------------------------------------------------------------------
    //
    // For each polycube vertex v we form an 8-bit "octant mask": bit
    // (dx + 2*dy + 4*dz) is 1 iff the cube with min-corner (v.x+dx-1,
    // v.y+dy-1, v.z+dz-1) is part of the polycube. (dx,dy,dz in {0,1}.)
    //
    // Up to the 48 cube symmetries (24 rotations + 24 reflections) the
    // possible non-empty / non-full configurations partition into the six
    // basic cases of paper Figure 5 (3, 4a, 4b, 5, 6a, 6b). We pick a
    // canonical mask for each and find the symmetry that brings the input
    // mask to that canonical one.
    //
    // Canonical octant masks (our choice; documented here):
    //   * Type 3  -- one of the 8 octants. Canonical: only (0,0,0) is in
    //                the polycube -> mask = 1 << 0 = 0x01.
    //                (This is "one corner cube" arrangement.)
    //   * Type 4a -- the four cubes lying on one side of the vertex (a slab).
    //                Canonical: dz=0 octants present -> bits 0,1,2,3 = 0x0F.
    //   * Type 4b -- two cubes sharing only an edge through v. Canonical:
    //                (0,0,0) and (1,1,0) -> bits 0 and 3 -> 0x09.
    //                (NOTE: this is the 4b "edge-pair" arrangement; under our
    //                dual partition the cell contains 4 facelets.)
    //   * Type 5  -- canonical: type-4a slab minus one corner. Bits 0,1,2 = 0x07
    //                (three coplanar cubes meeting at v).
    //
    // (Types 6a/6b are unsupported; we warn and treat them as empty.)

    private const byte MASK_TYPE3 = 0x01; // {(0,0,0)}
    private const byte MASK_TYPE4A = 0x0F; // {(0,0,0),(1,0,0),(0,1,0),(1,1,0)}
    private const byte MASK_TYPE4B = 0x09; // {(0,0,0),(1,1,0)}
    private const byte MASK_TYPE5 = 0x07;  // {(0,0,0),(1,0,0),(0,1,0)}

    public enum CellType { Empty = -1, Type3 = 0, Type4a = 1, Type4b = 2, Type5 = 3, Unsupported = 4 }

    // A signed-permutation matrix represented as 3 ints: each int is
    // axisIndex * 2 + (sign==-1 ? 1 : 0). axisIndex in {0,1,2}.
    public struct SymMat
    {
        public int Ax, Ay, Az; // column for output x,y,z
        public int Sx, Sy, Sz; // signs (+1 / -1)

        public Vector3 Apply(Vector3 v)
        {
            var arr = new[] { v.x, v.y, v.z };
            return new Vector3(Sx * arr[Ax], Sy * arr[Ay], Sz * arr[Az]);
        }

        // Apply to integer offset (dx,dy,dz) in {-1, 0, +1}.
        public Vector3Int ApplyInt(Vector3Int v)
        {
            var arr = new[] { v.x, v.y, v.z };
            return new Vector3Int(Sx * arr[Ax], Sy * arr[Ay], Sz * arr[Az]);
        }

        public override int GetHashCode()
        {
            return ((Ax * 2 + (Sx < 0 ? 1 : 0)) << 0)
                 | ((Ay * 2 + (Sy < 0 ? 1 : 0)) << 4)
                 | ((Az * 2 + (Sz < 0 ? 1 : 0)) << 8);
        }
        public override bool Equals(object o)
        {
            if (!(o is SymMat s)) return false;
            return Ax == s.Ax && Ay == s.Ay && Az == s.Az && Sx == s.Sx && Sy == s.Sy && Sz == s.Sz;
        }

        public static SymMat Identity => new SymMat { Ax = 0, Ay = 1, Az = 2, Sx = 1, Sy = 1, Sz = 1 };
    }

    private static List<SymMat> EnumerateAllSymmetries()
    {
        // 48 signed permutations of (x,y,z): 6 axis permutations * 8 sign choices.
        // Ordered with proper rotations (det = +1) first so that classification
        // prefers them; this matters because the shader's 5-op rotation table
        // only spans the 24 rotations -- reflections cannot be represented.
        var rots = new List<SymMat>(24);
        var refls = new List<SymMat>(24);
        var perms = new[]
        {
            new[]{0,1,2}, new[]{0,2,1}, new[]{1,0,2},
            new[]{1,2,0}, new[]{2,0,1}, new[]{2,1,0},
        };
        var permParity = new[] { 1, -1, -1, 1, 1, -1 };
        for (var p = 0; p < perms.Length; p++)
        {
            for (var s = 0; s < 8; s++)
            {
                var sx = (s & 1) != 0 ? -1 : 1;
                var sy = (s & 2) != 0 ? -1 : 1;
                var sz = (s & 4) != 0 ? -1 : 1;
                var det = permParity[p] * sx * sy * sz;
                var m = new SymMat { Ax = perms[p][0], Ay = perms[p][1], Az = perms[p][2], Sx = sx, Sy = sy, Sz = sz };
                if (det > 0) rots.Add(m); else refls.Add(m);
            }
        }
        rots.AddRange(refls);
        return rots;
    }

    private static readonly List<SymMat> AllSymmetries = EnumerateAllSymmetries();

    // Given an 8-bit octant mask, return the rotated mask after applying `m`.
    // Each bit (dx + 2*dy + 4*dz) corresponds to cube offset (dx, dy, dz) in
    // {0,1}^3 around the vertex. We map to centered offsets in {-1,+1}, apply
    // the symmetry, and remap back.
    private static byte RotateMask(byte mask, SymMat m)
    {
        byte outMask = 0;
        for (var bit = 0; bit < 8; bit++)
        {
            if ((mask & (1 << bit)) == 0) continue;
            var dx = (bit & 1);
            var dy = (bit >> 1) & 1;
            var dz = (bit >> 2) & 1;
            // Map {0,1} -> {-1,+1}
            var v = new Vector3Int(dx * 2 - 1, dy * 2 - 1, dz * 2 - 1);
            var r = m.ApplyInt(v);
            var ox = (r.x + 1) / 2;
            var oy = (r.y + 1) / 2;
            var oz = (r.z + 1) / 2;
            outMask |= (byte)(1 << (ox + 2 * oy + 4 * oz));
        }
        return outMask;
    }

    // Find a symmetry s such that RotateMask(inputMask, s) == canonicalMask
    // for some supported case. Returns (type, symmetry) or Empty/Unsupported.
    public static (CellType type, SymMat sym) ClassifyVertex(byte mask)
    {
        if (mask == 0 || mask == 0xFF) return (CellType.Empty, SymMat.Identity);

        var canonicals = new (CellType, byte)[]
        {
            (CellType.Type3,  MASK_TYPE3),
            (CellType.Type4a, MASK_TYPE4A),
            (CellType.Type4b, MASK_TYPE4B),
            (CellType.Type5,  MASK_TYPE5),
        };

        foreach (var (type, canonMask) in canonicals)
        {
            foreach (var s in AllSymmetries)
            {
                if (RotateMask(mask, s) == canonMask)
                    return (type, s);
            }
        }
        return (CellType.Unsupported, SymMat.Identity);
    }

    // ----------------------------------------------------------------------
    // The five rotation operations of paper §5.1, applied in shader order
    // (op1 first, then op2, op3, op4, op5). The shader uses the
    // post-op state of `sub` as the input to the next op.
    //
    //  op1: (r,s,t) -> (-t,-s,-r)
    //  op2: (r,s,t) -> ( s, t, r)
    //  op3: (r,s,t) -> ( s, t, r)
    //  op4: (r,s,t) -> (-r,-s, t)
    //  op5: (r,s,t) -> ( r,-s,-t)
    // ----------------------------------------------------------------------
    private static SymMat ApplyOp(SymMat m, int op)
    {
        // Compose m' = op_matrix * m so that (op_matrix * m)*v matches the
        // sequential transform on v. We apply each op to the *result columns*.
        // Easiest: apply op to the basis output (Ax,Sx),(Ay,Sy),(Az,Sz).
        Vector3Int X = AxisVec(m.Ax, m.Sx);
        Vector3Int Y = AxisVec(m.Ay, m.Sy);
        Vector3Int Z = AxisVec(m.Az, m.Sz);
        Vector3Int nx = X, ny = Y, nz = Z;
        switch (op)
        {
            case 1: // (r,s,t) -> (-t,-s,-r)
                nx = -Z; ny = -Y; nz = -X; break;
            case 2: // (r,s,t) -> (s,t,r)
            case 3:
                nx = Y; ny = Z; nz = X; break;
            case 4: // (r,s,t) -> (-r,-s,t)
                nx = -X; ny = -Y; nz = Z; break;
            case 5: // (r,s,t) -> (r,-s,-t)
                nx = X; ny = -Y; nz = -Z; break;
        }
        return FromAxes(nx, ny, nz);
    }

    private static Vector3Int AxisVec(int axis, int sign)
    {
        var v = Vector3Int.zero;
        if (axis == 0) v.x = sign;
        else if (axis == 1) v.y = sign;
        else v.z = sign;
        return v;
    }

    private static SymMat FromAxes(Vector3Int X, Vector3Int Y, Vector3Int Z)
    {
        // Each of X,Y,Z is one of +/- the three unit vectors.
        var m = new SymMat();
        (m.Ax, m.Sx) = AxisOf(X);
        (m.Ay, m.Sy) = AxisOf(Y);
        (m.Az, m.Sz) = AxisOf(Z);
        return m;
    }

    private static (int axis, int sign) AxisOf(Vector3Int v)
    {
        if (v.x != 0) return (0, v.x);
        if (v.y != 0) return (1, v.y);
        return (2, v.z);
    }

    // Build a dictionary mapping a fully-applied symmetry (the composition of
    // ops gated by a 5-bit pattern) to the 5-bit pattern that produces it.
    //
    // Bit numbering used here -- "rotBits":
    //   bit 0 -> op1 (gated by decoded.z, mask 32 in the V*32 frac test)
    //   bit 1 -> op2 (gated by decoded.w >= 0.234375 in frac(V*8))
    //   bit 2 -> op3 (gated by decoded.w >= 0.609375 in frac(V*8))
    //   bit 3 -> op4 (gated by decoded.x in frac(V*128))
    //   bit 4 -> op5 (gated by decoded.y in frac(V*64))
    //
    // (See EncodeByte for how those gating thresholds map into raw byte bits.)
    private static readonly Dictionary<SymMat, int> SymToRotBits = BuildSymTable();

    private static Dictionary<SymMat, int> BuildSymTable()
    {
        var d = new Dictionary<SymMat, int>();
        for (var bits = 0; bits < 32; bits++)
        {
            var m = SymMat.Identity;
            // Shader applies in order: op1, op2, op3, op4, op5.
            if ((bits & (1 << 0)) != 0) m = ApplyOp(m, 1);
            if ((bits & (1 << 1)) != 0) m = ApplyOp(m, 2);
            if ((bits & (1 << 2)) != 0) m = ApplyOp(m, 3);
            if ((bits & (1 << 3)) != 0) m = ApplyOp(m, 4);
            if ((bits & (1 << 4)) != 0) m = ApplyOp(m, 5);
            if (!d.ContainsKey(m)) d[m] = bits;
        }
        return d;
    }

    // The classification produces a symmetry that maps the *input* mask to
    // the *canonical* mask. The shader's rotation is supposed to map the
    // canonical orientation TO the input orientation (it rotates the
    // subcell-relative fragment position so that the projection formula
    // written for the canonical case applies). Hence we want the inverse of
    // the classifying symmetry. For sign-permutation matrices the inverse
    // is the transpose-with-signs.
    public static SymMat Inverse(SymMat m)
    {
        var inv = new SymMat();
        // Output axis Ax with sign Sx means input axis Ax goes to output
        // axis 0 (x) with that sign. So inverse: output axis n is input
        // axis whose Ai == n, with sign Si.
        for (var outAxis = 0; outAxis < 3; outAxis++)
        {
            int srcAxis = -1, sign = 1;
            if (m.Ax == outAxis) { srcAxis = 0; sign = m.Sx; }
            else if (m.Ay == outAxis) { srcAxis = 1; sign = m.Sy; }
            else if (m.Az == outAxis) { srcAxis = 2; sign = m.Sz; }
            if (outAxis == 0) { inv.Ax = srcAxis; inv.Sx = sign; }
            else if (outAxis == 1) { inv.Ay = srcAxis; inv.Sy = sign; }
            else { inv.Az = srcAxis; inv.Sz = sign; }
        }
        return inv;
    }

    // ----------------------------------------------------------------------
    // 3) Encode case + rotation as a single byte (the LUT's blue channel).
    // ----------------------------------------------------------------------
    //
    // The case is gated by the magnitude of V = b/255:
    //   V in [0,         31.55/255) -> case 3  -> byte [0,   31]
    //   V in [31.55/255, 63.43/255) -> case 4a -> byte [32,  63]
    //   V in [63.43/255, 95.30/255) -> case 4b -> byte [64,  95]
    //   V in [95.30/255, 127.18/255)-> case 5  -> byte [96, 127]
    //
    // The 5 rotation bits are extracted from V via the thresholds the shader
    // uses; we just brute-force the byte that matches a desired rotation
    // pattern within the case's byte range.

    public static int DecodeRotBits(byte b)
    {
        var V = b / 255f;
        var bits = 0;
        if (Frac(V * 128f) > 0.49f) bits |= (1 << 3); // op4
        if (Frac(V * 64f)  > 0.49f) bits |= (1 << 4); // op5
        if (Frac(V * 32f)  > 0.49f) bits |= (1 << 0); // op1
        var f8 = Frac(V * 8f);
        if (f8 > 0.234375f) bits |= (1 << 1);          // op2
        if (f8 > 0.609375f) bits |= (1 << 2);          // op3
        return bits;
    }

    public static CellType DecodeCase(byte b)
    {
        var V = b / 255f;
        if (V < 0.12375f) return CellType.Type3;
        if (V < 0.24875f) return CellType.Type4a;
        if (V < 0.37375f) return CellType.Type4b;
        if (V < 0.49875f) return CellType.Type5;
        return CellType.Unsupported;
    }

    private static float Frac(float x) { return x - Mathf.Floor(x); }

    public static byte EncodeByte(CellType type, int rotBits)
    {
        int lo, hi;
        switch (type)
        {
            case CellType.Type3:  lo = 0;  hi = 31;  break;
            case CellType.Type4a: lo = 32; hi = 63;  break;
            case CellType.Type4b: lo = 64; hi = 95;  break;
            case CellType.Type5:  lo = 96; hi = 127; break;
            default: throw new InvalidOperationException($"Unsupported cell type {type} for encoding.");
        }
        for (var b = lo; b <= hi; b++)
        {
            if (DecodeCase((byte)b) == type && DecodeRotBits((byte)b) == rotBits)
                return (byte)b;
        }
        throw new InvalidOperationException(
            $"No byte encoding exists for type {type} with rotBits=0x{rotBits:X2} " +
            "(no value in the case's byte range decodes to that rotation pattern).");
    }

    public static byte DecodeByte(byte b, out CellType type, out int rotBits)
    {
        type = DecodeCase(b);
        rotBits = DecodeRotBits(b);
        return b;
    }

    // ----------------------------------------------------------------------
    // 4) 3D-UV computation (paper §6.1). We perform the "warp + project"
    // step but skip the iterative MIPS optimization; the closest-surface-
    // point projection is a deliberate simplification of the paper's normal-
    // direction projection, chosen because it cannot fold over and works
    // robustly without optimization.
    // ----------------------------------------------------------------------

    public static Vector3[] ComputeTextureCoordinates(Mesh mesh, PolyCube polycube)
    {
        var meshBounds = mesh.bounds;
        var pcBounds = polycube.VertexBounds();

        // Build a similarity warp: uniform scale + translation that maps
        // mesh AABB center to polycube AABB center and shrinks to fit inside
        // the polycube AABB.
        var meshSize = meshBounds.size;
        var pcSize = pcBounds.size;
        var scale = Mathf.Min(
            pcSize.x / Mathf.Max(meshSize.x, 1e-8f),
            Mathf.Min(pcSize.y / Mathf.Max(meshSize.y, 1e-8f),
                      pcSize.z / Mathf.Max(meshSize.z, 1e-8f)));

        var meshVerts = mesh.vertices;
        var result = new Vector3[meshVerts.Length];
        for (var i = 0; i < meshVerts.Length; i++)
        {
            var warped = (meshVerts[i] - meshBounds.center) * scale + pcBounds.center;
            result[i] = ClosestPointOnPolycubeSurface(warped, polycube);
        }
        return result;
    }

    // Closest point on the polycube surface to p. Since the polycube is a
    // union of axis-aligned unit cubes, the surface is the union of the
    // outward-facing axis-aligned unit squares. We collect all *boundary*
    // faces (faces of a cube not shared with another cube in the polycube)
    // and find the closest point on any of them.
    private static Vector3 ClosestPointOnPolycubeSurface(Vector3 p, PolyCube polycube)
    {
        var best = Vector3.zero;
        var bestD2 = float.PositiveInfinity;

        foreach (var c in polycube.Cubes)
        {
            // 6 axis-aligned faces of this cube.
            for (var axis = 0; axis < 3; axis++)
            {
                for (var side = 0; side < 2; side++)
                {
                    var n = Vector3Int.zero;
                    if (axis == 0) n.x = side == 0 ? -1 : 1;
                    else if (axis == 1) n.y = side == 0 ? -1 : 1;
                    else n.z = side == 0 ? -1 : 1;

                    // Neighbor cube on the other side; if present, this face
                    // is interior and not part of the surface.
                    var neighbor = new Vector3Int(c.x + n.x, c.y + n.y, c.z + n.z);
                    if (polycube.ContainsCube(neighbor)) continue;

                    // Face rectangle: clamp p onto it.
                    var fMin = (Vector3)c;
                    var fMax = (Vector3)c + Vector3.one;
                    // Snap the axis coordinate to the face plane.
                    var planeCoord = (side == 0 ? c[axis] : c[axis] + 1);
                    var q = p;
                    q[axis] = planeCoord;
                    // Clamp other axes to face extents.
                    for (var a = 0; a < 3; a++)
                    {
                        if (a == axis) continue;
                        q[a] = Mathf.Clamp(p[a], fMin[a], fMax[a]);
                    }
                    var d2 = (q - p).sqrMagnitude;
                    if (d2 < bestD2) { bestD2 = d2; best = q; }
                }
            }
        }
        return best;
    }

    // TODO: Implement the iterative extended-MIPS optimization described in
    // paper §6.1 (last paragraph) -- per-vertex gradient descent on a
    // deformation energy with line search, followed by reprojection via P.
    // This is left as future work; the closest-point initialization above
    // is sufficient for the simple demo meshes shipped with this baker.
    public static void OptimizeParameterization(Mesh mesh, Vector3[] uvs, PolyCube polycube)
    {
        // Intentionally not implemented. See paper §6.1.
    }

    // ----------------------------------------------------------------------
    // 5) LUT + texture baking
    // ----------------------------------------------------------------------

    public class PatchInfo
    {
        public Vector3Int Cell;        // polycube vertex this entry corresponds to
        public CellType Type;
        public int RotBits;
        public byte EncodedByte;
        public int PatchPxX;           // top-left in pixels
        public int PatchPxY;
        public int PatchSquareletsW;   // size in squarelets (3)
        public int PatchSquareletsH;   // size in squarelets (2)
    }

    // We use a uniform 3x2-squarelet patch bounding rect for every supported
    // case. This is wasteful for type 3 (only 3 squarelets needed) and type 5
    // (5 needed) but keeps packing trivial and is documented as a deliberate
    // simplification. The exact squarelet-to-facelet mapping inside a patch is
    // not consumed by our demo content (we fill patches with a checkerboard),
    // so the per-cell color is what gets sampled in practice.
    private const int PATCH_SQUARELETS_W = 3;
    private const int PATCH_SQUARELETS_H = 2;

    public static (Texture2D tex, Dictionary<Vector3Int, PatchInfo> patches)
        BakeTexture(PolyCube polycube, int squareletSize, int textureWidth, int textureHeight)
    {
        var tex = new Texture2D(textureWidth, textureHeight, TextureFormat.RGBA32, false);
        tex.filterMode = FilterMode.Point;
        tex.wrapMode = TextureWrapMode.Clamp;

        // Fill with magenta so unused texels are visibly wrong.
        var clearPx = new Color32[textureWidth * textureHeight];
        var magenta = new Color32(255, 0, 255, 255);
        for (var i = 0; i < clearPx.Length; i++) clearPx[i] = magenta;
        tex.SetPixels32(clearPx);

        var patches = new Dictionary<Vector3Int, PatchInfo>();
        var verts = polycube.EnumerateVertices();

        // Compute LUT layout. LUT pixel (x,y) for cell (cx,cy,cz):
        //   x = cx + 16 * cz, y = cy   (matches shader's `dot(map, {1,0,16,0})`)
        // We need to reserve enough rows at the top for the LUT. Find max
        // cell coords first.
        var maxCellX = 0; var maxCellY = 0; var maxCellZ = 0;
        foreach (var v in verts)
        {
            if (v.x > maxCellX) maxCellX = v.x;
            if (v.y > maxCellY) maxCellY = v.y;
            if (v.z > maxCellZ) maxCellZ = v.z;
        }
        var lutMaxPxX = maxCellX + 16 * maxCellZ;
        var lutMaxPxY = maxCellY;
        if (lutMaxPxX >= textureWidth || lutMaxPxY >= textureHeight)
            throw new InvalidOperationException(
                $"LUT footprint ({lutMaxPxX + 1}x{lutMaxPxY + 1}) does not fit in texture " +
                $"({textureWidth}x{textureHeight}). Increase texture size.");

        // Patches go in a region below the LUT, starting at y = lutHeight
        // aligned up to squareletSize.
        var lutReservedRows = ((lutMaxPxY + 1) + squareletSize - 1) / squareletSize * squareletSize;

        // Pack patches left-to-right, top-to-bottom in squarelet units.
        var patchPxW = PATCH_SQUARELETS_W * squareletSize;
        var patchPxH = PATCH_SQUARELETS_H * squareletSize;
        var patchesPerRow = Mathf.Max(1, textureWidth / patchPxW);

        var assignedCount = 0;
        var typeCounts = new Dictionary<CellType, int>();

        // Stable iteration order: sort by (z,y,x).
        var sortedVerts = new List<Vector3Int>(verts);
        sortedVerts.Sort((a, b) =>
        {
            if (a.z != b.z) return a.z.CompareTo(b.z);
            if (a.y != b.y) return a.y.CompareTo(b.y);
            return a.x.CompareTo(b.x);
        });

        foreach (var v in sortedVerts)
        {
            // Build octant mask around v.
            byte mask = 0;
            for (var dz = 0; dz <= 1; dz++)
            for (var dy = 0; dy <= 1; dy++)
            for (var dx = 0; dx <= 1; dx++)
            {
                var c = new Vector3Int(v.x + dx - 1, v.y + dy - 1, v.z + dz - 1);
                if (polycube.ContainsCube(c)) mask |= (byte)(1 << (dx + 2 * dy + 4 * dz));
            }

            var (type, sym) = ClassifyVertex(mask);
            if (!typeCounts.ContainsKey(type)) typeCounts[type] = 0;
            typeCounts[type]++;

            if (type == CellType.Empty)
            {
                // Empty cell: leave LUT pixel as background; shader should not
                // sample here if 3D UVs stay on the surface.
                continue;
            }
            if (type == CellType.Unsupported)
            {
                Debug.LogWarning($"PolyCubeMapBaker: vertex {v} has unsupported configuration " +
                                 $"(mask=0x{mask:X2}, popcount={Popcount(mask)}). Treating as empty.");
                continue;
            }

            // The classifying symmetry rotates the actual cube arrangement
            // onto the canonical one. The same rotation, applied to the
            // fragment-relative position `sub`, brings the fragment into the
            // canonical frame where the projection formulas were derived --
            // which is exactly what the shader's e.R step does. So we want
            // `sym` itself (not its inverse) in the rotation table.
            var shaderSym = sym;
            if (!SymToRotBits.TryGetValue(shaderSym, out var rotBits))
            {
                // The 24 rotations are guaranteed reachable; reflections are
                // not. Polycube vertex classifications only need rotations
                // for the supported cases, but if a reflection sneaks in we
                // fall back to the identity-bits and warn.
                Debug.LogWarning($"PolyCubeMapBaker: symmetry for {v} not reachable by 5-bit ops " +
                                 "(likely a reflection). Falling back to 0.");
                rotBits = 0;
            }
            var encoded = EncodeByte(type, rotBits);

            // Patch placement.
            var slot = assignedCount;
            var patchCol = slot % patchesPerRow;
            var patchRow = slot / patchesPerRow;
            var patchPxX = patchCol * patchPxW;
            var patchPxY = lutReservedRows + patchRow * patchPxH;
            if (patchPxY + patchPxH > textureHeight)
                throw new InvalidOperationException(
                    "Ran out of texture space while packing patches; enlarge texture.");

            assignedCount++;

            // Write LUT pixel. R/G = patch position in squarelet units (so
            // that R*255*S/texWidth == patchPxX/texWidth -- see PolyCubeMap.cs).
            var lutX = v.x + 16 * v.z;
            var lutY = v.y;
            var lutR = (byte)(patchPxX / squareletSize);
            var lutG = (byte)(patchPxY / squareletSize);
            // Use byte-exact path so the LUT byte is preserved precisely.
            tex.SetPixels32(lutX, lutY, 1, 1, new[] { new Color32(lutR, lutG, encoded, 255) });

            // Fill the patch with a per-cell base color over a 2-color checker
            // squarelet pattern. Each squarelet is one solid color; alternating
            // squarelets darken slightly so the patch structure is visible.
            var baseCol = HsvHash(v);
            for (var sy = 0; sy < PATCH_SQUARELETS_H; sy++)
            {
                for (var sx = 0; sx < PATCH_SQUARELETS_W; sx++)
                {
                    var dark = ((sx + sy) & 1) == 1;
                    var col = dark ? new Color(baseCol.r * 0.6f, baseCol.g * 0.6f, baseCol.b * 0.6f, 1f) : baseCol;
                    var c32 = new Color32((byte)(col.r * 255), (byte)(col.g * 255), (byte)(col.b * 255), 255);
                    var px = patchPxX + sx * squareletSize;
                    var py = patchPxY + sy * squareletSize;
                    var block = new Color32[squareletSize * squareletSize];
                    for (var i = 0; i < block.Length; i++) block[i] = c32;
                    tex.SetPixels32(px, py, squareletSize, squareletSize, block);
                }
            }

            patches[v] = new PatchInfo
            {
                Cell = v,
                Type = type,
                RotBits = rotBits,
                EncodedByte = encoded,
                PatchPxX = patchPxX,
                PatchPxY = patchPxY,
                PatchSquareletsW = PATCH_SQUARELETS_W,
                PatchSquareletsH = PATCH_SQUARELETS_H,
            };
        }

        FillEmptyLUTEntries(
            (x, y, c) => tex.SetPixels32(x, y, 1, 1, new[] { c }),
            patches, maxCellX, maxCellY, maxCellZ, squareletSize);

        tex.Apply(false, false);

        // Brief utilization summary.
        Debug.Log($"PolyCubeMapBaker: baked {assignedCount} patch(es). " +
                  $"Type counts: " + FormatTypeCounts(typeCounts) +
                  $". Texture: {textureWidth}x{textureHeight}, squareletSize={squareletSize}, " +
                  $"patchPxSize={patchPxW}x{patchPxH}, patches/row={patchesPerRow}.");

        return (tex, patches);
    }

    private static int Popcount(byte b)
    {
        var n = 0; for (var i = 0; i < 8; i++) if ((b & (1 << i)) != 0) n++; return n;
    }

    // For voxelized / multi-cube polycubes, many cells in the LUT footprint
    // have no patch (interior vertices with mask=0xFF, empty exterior cells,
    // or unsupported 6a/6b configurations). Their LUT pixels would stay
    // magenta and fragments that land there in the shader produce visible
    // pink/garbage because the pink-as-bytes patch offset reads way off-atlas.
    //
    // Mitigation: copy the nearest non-empty cell's LUT entry into each
    // unwritten cell. The projection in that cell will then sample a nearby
    // patch -- not pixel-perfect, but no magenta. The proper fix (per paper
    // §3.2) is for the shader to project off-surface fragments back to T3,
    // which our runtime shader doesn't implement.
    private static void FillEmptyLUTEntries(
        System.Action<int, int, Color32> writePixel,
        Dictionary<Vector3Int, PatchInfo> patches,
        int maxCellX, int maxCellY, int maxCellZ,
        int squareletSize)
    {
        if (patches.Count == 0) return;
        var patchKeys = new List<Vector3Int>(patches.Keys);
        var filled = 0;
        for (var cz = 0; cz <= maxCellZ; cz++)
        for (var cy = 0; cy <= maxCellY; cy++)
        for (var cx = 0; cx <= maxCellX; cx++)
        {
            var key = new Vector3Int(cx, cy, cz);
            if (patches.ContainsKey(key)) continue;

            var bestD2 = int.MaxValue;
            PatchInfo bestPatch = null;
            foreach (var k in patchKeys)
            {
                var dx = k.x - cx; var dy = k.y - cy; var dz = k.z - cz;
                var d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < bestD2) { bestD2 = d2; bestPatch = patches[k]; }
            }
            if (bestPatch == null) continue;

            var lutX = cx + 16 * cz;
            var lutY = cy;
            var lutR = (byte)(bestPatch.PatchPxX / squareletSize);
            var lutG = (byte)(bestPatch.PatchPxY / squareletSize);
            writePixel(lutX, lutY, new Color32(lutR, lutG, bestPatch.EncodedByte, 255));
            filled++;
        }
        if (filled > 0)
            Debug.Log($"PolyCubeMapBaker: filled {filled} empty LUT cell(s) with nearest non-empty patch " +
                      "(off-surface fragments now sample an approximate neighbor instead of magenta).");
    }

    private static string FormatTypeCounts(Dictionary<CellType, int> d)
    {
        var sb = new System.Text.StringBuilder();
        foreach (var kv in d) { if (sb.Length > 0) sb.Append(", "); sb.Append($"{kv.Key}={kv.Value}"); }
        return sb.ToString();
    }

    // Deterministic, distinct-ish hue per cell.
    private static Color HsvHash(Vector3Int v)
    {
        unchecked
        {
            var h = (uint)(v.x * 73856093 ^ v.y * 19349663 ^ v.z * 83492791);
            var hue = (h & 0xFFFF) / 65535f;
            return Color.HSVToRGB(hue, 0.7f, 0.9f);
        }
    }

    // ----------------------------------------------------------------------
    // 6) End-to-end menu items
    // ----------------------------------------------------------------------

    private const string OutputDir = "Assets/Samples/Generated";
    private const int DefaultSquareletSize = 16;
    private const int DefaultTexWidth = 1024;
    private const int DefaultTexHeight = 1024;

    [MenuItem("PolyCubeMap/Bake/Sphere Sample (single-cube polycube)")]
    public static void BakeSphereSample()
    {
        var pc = new PolyCube(new[] { new Vector3Int(0, 0, 0) });
        var mesh = BuildIcosphereMesh(2);
        BakeAndPlace(pc, mesh, "sphere_pc");
    }

    [MenuItem("PolyCubeMap/Bake/L-shape Sample")]
    public static void BakeLShapeSample()
    {
        var pc = new PolyCube(new[]
        {
            new Vector3Int(0, 0, 0),
            new Vector3Int(1, 0, 0),
            new Vector3Int(1, 1, 0),
        });
        var mesh = BuildPolycubeSurfaceMesh(pc);
        BakeAndPlace(pc, mesh, "lshape_pc");
    }

    private static void BakeAndPlace(PolyCube pc, Mesh mesh, string baseName)
    {
        // 3D UVs (paper §6.1).
        var uvs3 = ComputeTextureCoordinates(mesh, pc);
        var uvList = new List<Vector3>(uvs3);
        mesh.SetUVs(0, uvList);
        mesh.UploadMeshData(false);

        var (tex, patches) = BakeTexture(pc, DefaultSquareletSize, DefaultTexWidth, DefaultTexHeight);

        // Ensure output directory exists.
        if (!Directory.Exists(OutputDir)) Directory.CreateDirectory(OutputDir);
        AssetDatabase.Refresh();

        // Save mesh asset.
        var meshPath = OutputDir + "/" + baseName + ".asset";
        var existingMesh = AssetDatabase.LoadAssetAtPath<Mesh>(meshPath);
        if (existingMesh != null) AssetDatabase.DeleteAsset(meshPath);
        AssetDatabase.CreateAsset(mesh, meshPath);

        // Save texture as PNG.
        var pngPath = OutputDir + "/" + baseName + "_texture.png";
        var bytes = tex.EncodeToPNG();
        File.WriteAllBytes(pngPath, bytes);
        AssetDatabase.ImportAsset(pngPath);

        // Re-import with Point filter & no compression so the LUT stays exact.
        var importer = (TextureImporter)AssetImporter.GetAtPath(pngPath);
        if (importer != null)
        {
            importer.filterMode = FilterMode.Point;
            importer.mipmapEnabled = false;
            importer.textureCompression = TextureImporterCompression.Uncompressed;
            importer.alphaIsTransparency = false;
            importer.npotScale = TextureImporterNPOTScale.None;
            importer.sRGBTexture = false;
            importer.SaveAndReimport();
        }
        AssetDatabase.SaveAssets();

        var savedMesh = AssetDatabase.LoadAssetAtPath<Mesh>(meshPath);
        var savedTex = AssetDatabase.LoadAssetAtPath<Texture2D>(pngPath);
        var mat = AssetDatabase.LoadAssetAtPath<Material>("Assets/Resources/PolyCubeMap.mat");

        // Create scene GameObject.
        var goName = baseName == "sphere_pc" ? "Sphere PolyCubeMap Sample" : "L-shape PolyCubeMap Sample";
        var existing = GameObject.Find(goName);
        if (existing != null) UnityEngine.Object.DestroyImmediate(existing);
        var go = new GameObject(goName);
        var mf = go.AddComponent<MeshFilter>();
        mf.sharedMesh = savedMesh;
        var mr = go.AddComponent<MeshRenderer>();
        if (mat != null) mr.sharedMaterial = mat;
        else Debug.LogWarning("PolyCubeMapBaker: PolyCubeMap.mat not found at Assets/Resources/PolyCubeMap.mat.");
        var pcm = go.AddComponent<PolyCubeMap>();
        pcm.PolyCubemapTexture = savedTex;
        pcm.SquareletSize = DefaultSquareletSize;

        Selection.activeGameObject = go;
        Debug.Log($"PolyCubeMapBaker: created '{goName}' in scene; mesh at {meshPath}, texture at {pngPath}, " +
                  $"{patches.Count} non-empty cells.");
    }

    // ----------------------------------------------------------------------
    // Simple procedural meshes for the demos.
    // ----------------------------------------------------------------------

    // Subdivided icosahedron, radius 1.
    private static Mesh BuildIcosphereMesh(int subdivisions)
    {
        // Initial icosahedron.
        var t = (1f + Mathf.Sqrt(5f)) * 0.5f;
        var verts = new List<Vector3>
        {
            new Vector3(-1,  t,  0), new Vector3( 1,  t,  0), new Vector3(-1, -t,  0), new Vector3( 1, -t,  0),
            new Vector3( 0, -1,  t), new Vector3( 0,  1,  t), new Vector3( 0, -1, -t), new Vector3( 0,  1, -t),
            new Vector3( t,  0, -1), new Vector3( t,  0,  1), new Vector3(-t,  0, -1), new Vector3(-t,  0,  1),
        };
        for (var i = 0; i < verts.Count; i++) verts[i] = verts[i].normalized;
        var tris = new List<int>
        {
            0,11,5, 0,5,1, 0,1,7, 0,7,10, 0,10,11,
            1,5,9, 5,11,4, 11,10,2, 10,7,6, 7,1,8,
            3,9,4, 3,4,2, 3,2,6, 3,6,8, 3,8,9,
            4,9,5, 2,4,11, 6,2,10, 8,6,7, 9,8,1,
        };

        var midpointCache = new Dictionary<long, int>();
        int Midpoint(int a, int b)
        {
            var key = a < b ? ((long)a << 32) | (uint)b : ((long)b << 32) | (uint)a;
            if (midpointCache.TryGetValue(key, out var idx)) return idx;
            var m = ((verts[a] + verts[b]) * 0.5f).normalized;
            verts.Add(m);
            idx = verts.Count - 1;
            midpointCache[key] = idx;
            return idx;
        }

        for (var s = 0; s < subdivisions; s++)
        {
            var nt = new List<int>(tris.Count * 4);
            for (var i = 0; i < tris.Count; i += 3)
            {
                var a = tris[i]; var b = tris[i + 1]; var c = tris[i + 2];
                var ab = Midpoint(a, b);
                var bc = Midpoint(b, c);
                var ca = Midpoint(c, a);
                nt.Add(a); nt.Add(ab); nt.Add(ca);
                nt.Add(b); nt.Add(bc); nt.Add(ab);
                nt.Add(c); nt.Add(ca); nt.Add(bc);
                nt.Add(ab); nt.Add(bc); nt.Add(ca);
            }
            tris = nt;
        }

        var mesh = new Mesh();
        mesh.vertices = verts.ToArray();
        mesh.triangles = tris.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
        return mesh;
    }

    // The polycube's own surface as a triangle mesh (boundary faces only).
    private static Mesh BuildPolycubeSurfaceMesh(PolyCube pc)
    {
        var verts = new List<Vector3>();
        var tris = new List<int>();

        foreach (var c in pc.Cubes)
        {
            for (var axis = 0; axis < 3; axis++)
            {
                for (var side = 0; side < 2; side++)
                {
                    var n = Vector3Int.zero;
                    if (axis == 0) n.x = side == 0 ? -1 : 1;
                    else if (axis == 1) n.y = side == 0 ? -1 : 1;
                    else n.z = side == 0 ? -1 : 1;
                    var neighbor = new Vector3Int(c.x + n.x, c.y + n.y, c.z + n.z);
                    if (pc.ContainsCube(neighbor)) continue;

                    // Emit a quad as two triangles.
                    var planeCoord = side == 0 ? c[axis] : c[axis] + 1;
                    var u = (axis + 1) % 3;
                    var v = (axis + 2) % 3;
                    Vector3 P(int au, int av)
                    {
                        var p = new Vector3();
                        p[axis] = planeCoord;
                        p[u] = c[u] + au;
                        p[v] = c[v] + av;
                        return p;
                    }
                    var i0 = verts.Count;
                    verts.Add(P(0, 0)); verts.Add(P(1, 0)); verts.Add(P(1, 1)); verts.Add(P(0, 1));
                    // Winding so normal points outward (along n).
                    if (side == 0)
                    {
                        tris.Add(i0 + 0); tris.Add(i0 + 2); tris.Add(i0 + 1);
                        tris.Add(i0 + 0); tris.Add(i0 + 3); tris.Add(i0 + 2);
                    }
                    else
                    {
                        tris.Add(i0 + 0); tris.Add(i0 + 1); tris.Add(i0 + 2);
                        tris.Add(i0 + 0); tris.Add(i0 + 2); tris.Add(i0 + 3);
                    }
                }
            }
        }

        var mesh = new Mesh();
        mesh.vertices = verts.ToArray();
        mesh.triangles = tris.ToArray();
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
        return mesh;
    }

    // ----------------------------------------------------------------------
    // 7) Self tests (run at static init). On failure, log an error but do
    // not throw -- the user can still open the editor.
    // ----------------------------------------------------------------------

    static PolyCubeMapBaker()
    {
        try
        {
            RunSelfTests();
        }
        catch (Exception e)
        {
            Debug.LogError($"PolyCubeMapBaker self-tests threw: {e}");
        }
    }

    // ----------------------------------------------------------------------
    // 8) Voxelization + arbitrary-mesh baking (new).
    //
    // These additions let a user pick any Unity mesh + an optional Texture2D
    // and bake a PolyCubeMap from it. The voxelization is a deliberately
    // simple triangle-sampling + flood-fill scheme that is good enough for
    // small demo meshes but is NOT what you'd ship: the production-quality
    // construction is the Fu/Bai/Liu 2016 polycube algorithm, see
    // Material/paper.pdf and the C++ in Material/*.cpp.
    // ----------------------------------------------------------------------

    // Shader constraint: LUT pixel is `cellX + 16 * cellZ`, so cellX must be
    // < 16. Vertex coords go up to resolution+1, so max safe resolution is 14.
    private const int MaxResolution = 14;
    private const int MinResolution = 2;

    public static PolyCube Voxelize(Mesh mesh, int resolution, bool fillInterior)
    {
        if (resolution > MaxResolution)
        {
            Debug.LogError($"PolyCubeMapBaker.Voxelize: resolution {resolution} exceeds shader cap " +
                           $"({MaxResolution}); clamping. The shader's `cellX + 16*cellZ` LUT packing " +
                           "collides for cellX >= 16.");
            resolution = MaxResolution;
        }
        if (resolution < MinResolution)
        {
            Debug.LogError($"PolyCubeMapBaker.Voxelize: resolution {resolution} below minimum " +
                           $"({MinResolution}); clamping.");
            resolution = MinResolution;
        }

        Vector3[] verts;
        int[] tris;
        try
        {
            verts = mesh.vertices;
            tris = mesh.triangles;
        }
        catch (Exception e)
        {
            Debug.LogError($"PolyCubeMapBaker.Voxelize: cannot read mesh '{mesh.name}' " +
                           "(is Read/Write enabled in the importer?). " + e.Message);
            return new PolyCube();
        }
        if (tris.Length == 0)
        {
            Debug.LogError("PolyCubeMapBaker.Voxelize: mesh has no triangles.");
            return new PolyCube();
        }

        var b = mesh.bounds;
        var size = b.size;
        var longest = Mathf.Max(size.x, Mathf.Max(size.y, size.z));
        if (longest <= 1e-8f)
        {
            Debug.LogError("PolyCubeMapBaker.Voxelize: mesh has zero extent.");
            return new PolyCube();
        }
        var voxelSize = longest / resolution;

        // Whole grid extents per axis (not padded to a cube).
        var nx = Mathf.Max(1, Mathf.CeilToInt(size.x / voxelSize));
        var ny = Mathf.Max(1, Mathf.CeilToInt(size.y / voxelSize));
        var nz = Mathf.Max(1, Mathf.CeilToInt(size.z / voxelSize));
        // The longest axis should land exactly on `resolution`; the rounding
        // above guards against tiny FP slop.
        nx = Mathf.Min(nx, MaxResolution);
        ny = Mathf.Min(ny, MaxResolution);
        nz = Mathf.Min(nz, MaxResolution);

        var origin = b.min;
        var shell = new bool[nx, ny, nz];

        Vector3Int ToVoxel(Vector3 p)
        {
            var local = (p - origin) / voxelSize;
            return new Vector3Int(
                Mathf.Clamp(Mathf.FloorToInt(local.x), 0, nx - 1),
                Mathf.Clamp(Mathf.FloorToInt(local.y), 0, ny - 1),
                Mathf.Clamp(Mathf.FloorToInt(local.z), 0, nz - 1));
        }

        // Triangle sampling: pick a step count from the longest edge so we
        // hit every voxel the triangle crosses (factor 0.4 < 0.5 to over-
        // sample slightly and avoid hairline gaps).
        for (var i = 0; i < tris.Length; i += 3)
        {
            var v0 = verts[tris[i]];
            var v1 = verts[tris[i + 1]];
            var v2 = verts[tris[i + 2]];
            var e0 = (v1 - v0).magnitude;
            var e1 = (v2 - v1).magnitude;
            var e2 = (v0 - v2).magnitude;
            var maxEdge = Mathf.Max(e0, Mathf.Max(e1, e2));
            var steps = Mathf.Max(1, Mathf.CeilToInt(maxEdge / (voxelSize * 0.4f)));
            for (var su = 0; su <= steps; su++)
            {
                var u = (float)su / steps;
                for (var sv = 0; sv <= steps - su; sv++)
                {
                    var w = (float)sv / steps;
                    var p = v0 + (v1 - v0) * u + (v2 - v0) * w;
                    var vox = ToVoxel(p);
                    shell[vox.x, vox.y, vox.z] = true;
                }
            }
        }

        var cubes = new List<Vector3Int>();

        if (fillInterior)
        {
            // 3D BFS flood from (-1,-1,-1) over the bbox padded by 1 in every
            // direction. A voxel inside the bbox that is neither shell nor
            // reached by the flood is interior. Will leak on open meshes --
            // we cannot detect that cheaply, so we just log a hint when the
            // interior count looks suspicious.
            var pnx = nx + 2; var pny = ny + 2; var pnz = nz + 2;
            var visited = new bool[pnx, pny, pnz];
            // Use packed int to avoid Vector3Int alloc pressure in the queue.
            var q = new Queue<int>();
            int Pack(int x, int y, int z) => (x * pny + y) * pnz + z;
            visited[0, 0, 0] = true;
            q.Enqueue(Pack(0, 0, 0));
            var dirs = new[]
            {
                new Vector3Int(1, 0, 0), new Vector3Int(-1, 0, 0),
                new Vector3Int(0, 1, 0), new Vector3Int(0, -1, 0),
                new Vector3Int(0, 0, 1), new Vector3Int(0, 0, -1),
            };
            while (q.Count > 0)
            {
                var packed = q.Dequeue();
                var z = packed % pnz;
                var y = (packed / pnz) % pny;
                var x = packed / (pnz * pny);
                foreach (var d in dirs)
                {
                    var px = x + d.x; var py = y + d.y; var pz = z + d.z;
                    if (px < 0 || py < 0 || pz < 0 || px >= pnx || py >= pny || pz >= pnz) continue;
                    if (visited[px, py, pz]) continue;
                    // Bbox-local coords (subtract the 1-pad).
                    var lx = px - 1; var ly = py - 1; var lz = pz - 1;
                    var inBbox = lx >= 0 && ly >= 0 && lz >= 0 && lx < nx && ly < ny && lz < nz;
                    if (inBbox && shell[lx, ly, lz]) continue; // blocked by shell
                    visited[px, py, pz] = true;
                    q.Enqueue(Pack(px, py, pz));
                }
            }
            for (var x = 0; x < nx; x++)
            for (var y = 0; y < ny; y++)
            for (var z = 0; z < nz; z++)
            {
                if (shell[x, y, z] || !visited[x + 1, y + 1, z + 1])
                    cubes.Add(new Vector3Int(x, y, z));
            }
            // Heuristic leak hint: if every voxel ended up filled, the flood
            // probably leaked through an open mesh hole.
            if (cubes.Count == nx * ny * nz)
                Debug.LogWarning("PolyCubeMapBaker.Voxelize: fillInterior produced a fully-solid grid; " +
                                 "the mesh may not be watertight (interior flood leaked).");
        }
        else
        {
            for (var x = 0; x < nx; x++)
            for (var y = 0; y < ny; y++)
            for (var z = 0; z < nz; z++)
                if (shell[x, y, z]) cubes.Add(new Vector3Int(x, y, z));
        }

        if (cubes.Count == 0)
        {
            Debug.LogWarning("PolyCubeMapBaker.Voxelize: no cubes produced; mesh too thin for grid?");
        }

        return new PolyCube(cubes);
    }

    public static (Texture2D tex, Dictionary<Vector3Int, PatchInfo> patches)
        BakeTextureFromMesh(PolyCube polycube, Mesh mesh, Vector3[] uvs3, Texture2D sourceTexture,
                            int squareletSize, int textureWidth, int textureHeight)
    {
        // If no source, fall back to the existing checkerboard baker.
        if (sourceTexture == null)
            return BakeTexture(polycube, squareletSize, textureWidth, textureHeight);

        // Attempt to read the source texture; if Read/Write is off, the call
        // throws and we fall back to the checkerboard with a warning.
        Color32[] srcCol = null;
        int sw = 0, sh = 0;
        try
        {
            srcCol = sourceTexture.GetPixels32();
            sw = sourceTexture.width;
            sh = sourceTexture.height;
        }
        catch (Exception e)
        {
            Debug.LogWarning($"PolyCubeMapBaker.BakeTextureFromMesh: cannot read source texture " +
                             $"'{sourceTexture.name}' (Read/Write disabled?): {e.Message}. " +
                             "Falling back to checkerboard fill.");
            return BakeTexture(polycube, squareletSize, textureWidth, textureHeight);
        }

        var meshUV = new List<Vector2>();
        mesh.GetUVs(0, meshUV);
        if (meshUV.Count != mesh.vertexCount)
        {
            Debug.LogWarning($"PolyCubeMapBaker.BakeTextureFromMesh: mesh '{mesh.name}' has " +
                             $"{meshUV.Count} UVs vs {mesh.vertexCount} vertices; falling back to " +
                             "checkerboard.");
            return BakeTexture(polycube, squareletSize, textureWidth, textureHeight);
        }

        var meshTris = mesh.triangles;
        var meshVerts = mesh.vertices;
        if (meshTris.Length == 0)
        {
            Debug.LogWarning("PolyCubeMapBaker.BakeTextureFromMesh: mesh has no triangles; " +
                             "falling back to checkerboard.");
            return BakeTexture(polycube, squareletSize, textureWidth, textureHeight);
        }

        // For huge meshes, skip the per-pixel triangle scan and fall back to
        // a vertex-only nearest search (much faster, slightly worse quality).
        var triCount = meshTris.Length / 3;
        var vertexOnlySampling = triCount > 20000;
        if (vertexOnlySampling)
            Debug.LogWarning($"PolyCubeMapBaker.BakeTextureFromMesh: mesh has {triCount} triangles " +
                             "(>20000); using vertex-only nearest sampling for speed.");

        // Build inverse warp: polycube-space p -> mesh-space q.
        // Forward: warped = (q - meshCenter) * scale + pcCenter
        // Inverse: q = (warped - pcCenter) / scale + meshCenter
        var meshBounds = mesh.bounds;
        var pcBounds = polycube.VertexBounds();
        var meshSize = meshBounds.size;
        var pcSize = pcBounds.size;
        var scale = Mathf.Min(
            pcSize.x / Mathf.Max(meshSize.x, 1e-8f),
            Mathf.Min(pcSize.y / Mathf.Max(meshSize.y, 1e-8f),
                      pcSize.z / Mathf.Max(meshSize.z, 1e-8f)));
        var invScale = 1f / Mathf.Max(scale, 1e-12f);

        Vector3 PolycubeToMesh(Vector3 p) => (p - pcBounds.center) * invScale + meshBounds.center;

        // Set up the atlas.
        var tex = new Texture2D(textureWidth, textureHeight, TextureFormat.RGBA32, false);
        tex.filterMode = FilterMode.Point;
        tex.wrapMode = TextureWrapMode.Clamp;
        var atlas = new Color32[textureWidth * textureHeight];
        var magenta = new Color32(255, 0, 255, 255);
        for (var i = 0; i < atlas.Length; i++) atlas[i] = magenta;

        var patches = new Dictionary<Vector3Int, PatchInfo>();
        var pcVerts = polycube.EnumerateVertices();
        var maxCellX = 0; var maxCellY = 0; var maxCellZ = 0;
        foreach (var v in pcVerts)
        {
            if (v.x > maxCellX) maxCellX = v.x;
            if (v.y > maxCellY) maxCellY = v.y;
            if (v.z > maxCellZ) maxCellZ = v.z;
        }
        var lutMaxPxX = maxCellX + 16 * maxCellZ;
        var lutMaxPxY = maxCellY;
        if (lutMaxPxX >= textureWidth || lutMaxPxY >= textureHeight)
            throw new InvalidOperationException(
                $"LUT footprint ({lutMaxPxX + 1}x{lutMaxPxY + 1}) does not fit in texture " +
                $"({textureWidth}x{textureHeight}). Increase texture size.");
        var lutReservedRows = ((lutMaxPxY + 1) + squareletSize - 1) / squareletSize * squareletSize;

        var patchPxW = PATCH_SQUARELETS_W * squareletSize;
        var patchPxH = PATCH_SQUARELETS_H * squareletSize;
        var patchesPerRow = Mathf.Max(1, textureWidth / patchPxW);

        var sortedVerts = new List<Vector3Int>(pcVerts);
        sortedVerts.Sort((a, b) =>
        {
            if (a.z != b.z) return a.z.CompareTo(b.z);
            if (a.y != b.y) return a.y.CompareTo(b.y);
            return a.x.CompareTo(b.x);
        });

        var assignedCount = 0;
        var typeCounts = new Dictionary<CellType, int>();

        foreach (var v in sortedVerts)
        {
            byte mask = 0;
            for (var dz = 0; dz <= 1; dz++)
            for (var dy = 0; dy <= 1; dy++)
            for (var dx = 0; dx <= 1; dx++)
            {
                var c = new Vector3Int(v.x + dx - 1, v.y + dy - 1, v.z + dz - 1);
                if (polycube.ContainsCube(c)) mask |= (byte)(1 << (dx + 2 * dy + 4 * dz));
            }

            var (type, sym) = ClassifyVertex(mask);
            if (!typeCounts.ContainsKey(type)) typeCounts[type] = 0;
            typeCounts[type]++;

            if (type == CellType.Empty) continue;
            if (type == CellType.Unsupported)
            {
                Debug.LogWarning($"PolyCubeMapBaker: vertex {v} has unsupported configuration " +
                                 $"(mask=0x{mask:X2}). Treating as empty.");
                continue;
            }

            if (!SymToRotBits.TryGetValue(sym, out var rotBits))
            {
                Debug.LogWarning($"PolyCubeMapBaker: symmetry for {v} not reachable by 5-bit ops.");
                rotBits = 0;
            }
            var encoded = EncodeByte(type, rotBits);
            var invSym = Inverse(sym);

            var slot = assignedCount;
            var patchCol = slot % patchesPerRow;
            var patchRow = slot / patchesPerRow;
            var patchPxX = patchCol * patchPxW;
            var patchPxY = lutReservedRows + patchRow * patchPxH;
            if (patchPxY + patchPxH > textureHeight)
                throw new InvalidOperationException(
                    "Ran out of texture space while packing patches; enlarge texture.");
            assignedCount++;

            // LUT byte: write directly into our backing atlas array.
            var lutX = v.x + 16 * v.z;
            var lutY = v.y;
            var lutR = (byte)(patchPxX / squareletSize);
            var lutG = (byte)(patchPxY / squareletSize);
            atlas[lutY * textureWidth + lutX] = new Color32(lutR, lutG, encoded, 255);

            // Fill each of the 6 squarelets. Uniform projection scheme:
            // squarelet (sx, sy) in [0..2]x[0..1] covers polycube-space offset
            //   ((sx-1) + px/(S-1)) * 0.5  on X
            //   ((sy-1) + py/(S-1)) * 0.5  on Y
            //   0                          on Z
            // from the cell center, then rotated by inv(sym) into world frame.
            // This is a deliberate simplification -- the shader's per-case
            // squarelet-to-facelet mapping is more nuanced, but for the
            // supported types the facelets the shader actually samples land
            // in roughly the right spots. The unused squarelets will be wrong
            // and won't be sampled.
            var center = (Vector3)v;
            for (var sy = 0; sy < PATCH_SQUARELETS_H; sy++)
            for (var sx = 0; sx < PATCH_SQUARELETS_W; sx++)
            {
                for (var py = 0; py < squareletSize; py++)
                {
                    for (var px = 0; px < squareletSize; px++)
                    {
                        // Map pixel within squarelet to local 2D offset.
                        var fx = squareletSize > 1 ? (float)px / (squareletSize - 1) : 0.5f;
                        var fy = squareletSize > 1 ? (float)py / (squareletSize - 1) : 0.5f;
                        var local = new Vector3(((sx - 1) + fx) * 0.5f,
                                                ((sy - 1) + fy) * 0.5f,
                                                0f);
                        var pPoly = center + invSym.Apply(local);
                        var pMesh = PolycubeToMesh(pPoly);

                        Vector2 uv;
                        if (vertexOnlySampling)
                            uv = NearestVertexUV(pMesh, meshVerts, meshUV);
                        else
                            uv = NearestTriangleUV(pMesh, meshVerts, meshTris, meshUV);

                        var su = Mathf.Clamp01(uv.x);
                        var sv = Mathf.Clamp01(uv.y);
                        var ix = Mathf.Clamp(Mathf.FloorToInt(su * sw), 0, sw - 1);
                        var iy = Mathf.Clamp(Mathf.FloorToInt(sv * sh), 0, sh - 1);
                        var col = srcCol[iy * sw + ix];

                        var atlasX = patchPxX + sx * squareletSize + px;
                        var atlasY = patchPxY + sy * squareletSize + py;
                        atlas[atlasY * textureWidth + atlasX] = col;
                    }
                }
            }

            patches[v] = new PatchInfo
            {
                Cell = v,
                Type = type,
                RotBits = rotBits,
                EncodedByte = encoded,
                PatchPxX = patchPxX,
                PatchPxY = patchPxY,
                PatchSquareletsW = PATCH_SQUARELETS_W,
                PatchSquareletsH = PATCH_SQUARELETS_H,
            };
        }

        FillEmptyLUTEntries(
            (x, y, c) => atlas[y * textureWidth + x] = c,
            patches, maxCellX, maxCellY, maxCellZ, squareletSize);

        tex.SetPixels32(atlas);
        tex.Apply(false, false);

        Debug.Log($"PolyCubeMapBaker: baked-from-mesh {assignedCount} patch(es). " +
                  $"Type counts: " + FormatTypeCounts(typeCounts) +
                  $". Texture: {textureWidth}x{textureHeight}, squareletSize={squareletSize}, " +
                  $"patchPxSize={patchPxW}x{patchPxH}, patches/row={patchesPerRow}.");

        return (tex, patches);
    }

    // Nearest-vertex fallback (used when the mesh is too big for per-triangle).
    private static Vector2 NearestVertexUV(Vector3 p, Vector3[] verts, List<Vector2> uvs)
    {
        var bestD2 = float.PositiveInfinity;
        var bestUV = Vector2.zero;
        for (var i = 0; i < verts.Length; i++)
        {
            var d2 = (verts[i] - p).sqrMagnitude;
            if (d2 < bestD2) { bestD2 = d2; bestUV = uvs[i]; }
        }
        return bestUV;
    }

    // Closest point on the triangle soup, interpolated UV at that point.
    private static Vector2 NearestTriangleUV(Vector3 p, Vector3[] verts, int[] tris, List<Vector2> uvs)
    {
        var bestD2 = float.PositiveInfinity;
        var bestUV = Vector2.zero;
        for (var t = 0; t < tris.Length; t += 3)
        {
            var ia = tris[t]; var ib = tris[t + 1]; var ic = tris[t + 2];
            var a = verts[ia]; var b = verts[ib]; var c = verts[ic];
            ClosestPointOnTriangle(p, a, b, c, out var q, out var u, out var v, out var w);
            var d2 = (q - p).sqrMagnitude;
            if (d2 < bestD2)
            {
                bestD2 = d2;
                bestUV = uvs[ia] * u + uvs[ib] * v + uvs[ic] * w;
            }
        }
        return bestUV;
    }

    // Closest-point-on-triangle (Ericson "Real-Time Collision Detection", §5.1.5).
    // Returns the point + barycentric coordinates (u for a, v for b, w for c).
    private static void ClosestPointOnTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c,
                                               out Vector3 q, out float u, out float v, out float w)
    {
        var ab = b - a;
        var ac = c - a;
        var ap = p - a;
        var d1 = Vector3.Dot(ab, ap);
        var d2 = Vector3.Dot(ac, ap);
        if (d1 <= 0f && d2 <= 0f) { q = a; u = 1f; v = 0f; w = 0f; return; }

        var bp = p - b;
        var d3 = Vector3.Dot(ab, bp);
        var d4 = Vector3.Dot(ac, bp);
        if (d3 >= 0f && d4 <= d3) { q = b; u = 0f; v = 1f; w = 0f; return; }

        var vc = d1 * d4 - d3 * d2;
        if (vc <= 0f && d1 >= 0f && d3 <= 0f)
        {
            var t = d1 / (d1 - d3);
            q = a + t * ab; u = 1f - t; v = t; w = 0f; return;
        }

        var cp = p - c;
        var d5 = Vector3.Dot(ab, cp);
        var d6 = Vector3.Dot(ac, cp);
        if (d6 >= 0f && d5 <= d6) { q = c; u = 0f; v = 0f; w = 1f; return; }

        var vb = d5 * d2 - d1 * d6;
        if (vb <= 0f && d2 >= 0f && d6 <= 0f)
        {
            var t = d2 / (d2 - d6);
            q = a + t * ac; u = 1f - t; v = 0f; w = t; return;
        }

        var va = d3 * d6 - d5 * d4;
        if (va <= 0f && (d4 - d3) >= 0f && (d5 - d6) >= 0f)
        {
            var t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            q = b + t * (c - b); u = 0f; v = 1f - t; w = t; return;
        }

        var denom = 1f / (va + vb + vc);
        v = vb * denom;
        w = vc * denom;
        u = 1f - v - w;
        q = a + ab * v + ac * w;
    }

    // Full pipeline used by the menu / window.
    public static GameObject BakeFromMesh(Mesh mesh, Texture2D sourceTexture, int resolution,
                                          bool fillInterior, int squareletSize,
                                          int textureWidth, int textureHeight, string outputBaseName)
    {
        if (mesh == null)
        {
            Debug.LogError("PolyCubeMapBaker.BakeFromMesh: mesh is null.");
            return null;
        }
        if (string.IsNullOrWhiteSpace(outputBaseName)) outputBaseName = mesh.name + "_pcm";

        var pc = Voxelize(mesh, resolution, fillInterior);
        if (pc.Cubes.Count == 0)
        {
            Debug.LogError("PolyCubeMapBaker.BakeFromMesh: voxelization produced no cubes; aborting.");
            return null;
        }

        // Compute 3D UVs from the polycube and a copy of the mesh (so we
        // don't stomp the user's source mesh asset).
        var bakedMesh = UnityEngine.Object.Instantiate(mesh);
        bakedMesh.name = outputBaseName;
        var uvs3 = ComputeTextureCoordinates(bakedMesh, pc);
        bakedMesh.SetUVs(0, new List<Vector3>(uvs3));
        bakedMesh.UploadMeshData(false);

        var (tex, patches) = sourceTexture != null
            ? BakeTextureFromMesh(pc, mesh, uvs3, sourceTexture, squareletSize, textureWidth, textureHeight)
            : BakeTexture(pc, squareletSize, textureWidth, textureHeight);

        if (!Directory.Exists(OutputDir)) Directory.CreateDirectory(OutputDir);
        AssetDatabase.Refresh();

        var meshPath = OutputDir + "/" + outputBaseName + ".asset";
        var existingMesh = AssetDatabase.LoadAssetAtPath<Mesh>(meshPath);
        if (existingMesh != null) AssetDatabase.DeleteAsset(meshPath);
        AssetDatabase.CreateAsset(bakedMesh, meshPath);

        var pngPath = OutputDir + "/" + outputBaseName + "_texture.png";
        var bytes = tex.EncodeToPNG();
        File.WriteAllBytes(pngPath, bytes);
        AssetDatabase.ImportAsset(pngPath);

        var importer = (TextureImporter)AssetImporter.GetAtPath(pngPath);
        if (importer != null)
        {
            importer.filterMode = FilterMode.Point;
            importer.mipmapEnabled = false;
            importer.textureCompression = TextureImporterCompression.Uncompressed;
            importer.alphaIsTransparency = false;
            importer.npotScale = TextureImporterNPOTScale.None;
            importer.sRGBTexture = false;
            importer.SaveAndReimport();
        }
        AssetDatabase.SaveAssets();

        var savedMesh = AssetDatabase.LoadAssetAtPath<Mesh>(meshPath);
        var savedTex = AssetDatabase.LoadAssetAtPath<Texture2D>(pngPath);
        var mat = AssetDatabase.LoadAssetAtPath<Material>("Assets/Resources/PolyCubeMap.mat");

        var goName = outputBaseName;
        var existing = GameObject.Find(goName);
        if (existing != null) UnityEngine.Object.DestroyImmediate(existing);
        var go = new GameObject(goName);
        var mf = go.AddComponent<MeshFilter>();
        mf.sharedMesh = savedMesh;
        var mr = go.AddComponent<MeshRenderer>();
        if (mat != null) mr.sharedMaterial = mat;
        else Debug.LogWarning("PolyCubeMapBaker: PolyCubeMap.mat not found at Assets/Resources/PolyCubeMap.mat.");
        var pcm = go.AddComponent<PolyCubeMap>();
        pcm.PolyCubemapTexture = savedTex;
        pcm.SquareletSize = squareletSize;

        Selection.activeGameObject = go;
        Debug.Log($"PolyCubeMapBaker.BakeFromMesh: '{goName}' baked from mesh '{mesh.name}' " +
                  $"({pc.Cubes.Count} cubes, {patches.Count} patches, res={resolution}, " +
                  $"fillInterior={fillInterior}, src={(sourceTexture != null ? sourceTexture.name : "none")}).");
        return go;
    }

    private static void RunSelfTests()
    {
        // (a) Byte roundtrip: decoding then re-encoding (in the same case)
        // should give back the same byte.
        for (var b = 0; b < 256; b++)
        {
            var type = DecodeCase((byte)b);
            if (type == CellType.Unsupported) continue; // case 6 region; ignore
            var rot = DecodeRotBits((byte)b);
            byte re;
            try { re = EncodeByte(type, rot); }
            catch (InvalidOperationException ex)
            {
                Debug.LogError($"PolyCubeMapBaker self-test: EncodeByte failed for b={b}: {ex.Message}");
                return;
            }
            // We only require re-decoded byte matches in (type, rot); the
            // exact byte may differ because multiple bytes can decode to the
            // same (type, rot) pair.
            if (DecodeCase(re) != type || DecodeRotBits(re) != rot)
            {
                Debug.LogError($"PolyCubeMapBaker self-test: roundtrip mismatch at b={b}: " +
                               $"got type={DecodeCase(re)}, rot={DecodeRotBits(re):X2} " +
                               $"(expected type={type}, rot={rot:X2}).");
                return;
            }
        }

        // (b) Single-cube polycube: 8 vertices, all type 3, all encodable.
        var pc = new PolyCube(new[] { new Vector3Int(0, 0, 0) });
        var verts = pc.EnumerateVertices();
        if (verts.Count != 8)
        {
            Debug.LogError($"PolyCubeMapBaker self-test: expected 8 vertices, got {verts.Count}.");
            return;
        }
        var t3Count = 0;
        foreach (var v in verts)
        {
            byte mask = 0;
            for (var dz = 0; dz <= 1; dz++)
            for (var dy = 0; dy <= 1; dy++)
            for (var dx = 0; dx <= 1; dx++)
            {
                var c = new Vector3Int(v.x + dx - 1, v.y + dy - 1, v.z + dz - 1);
                if (pc.ContainsCube(c)) mask |= (byte)(1 << (dx + 2 * dy + 4 * dz));
            }
            var (type, sym) = ClassifyVertex(mask);
            if (type != CellType.Type3)
            {
                Debug.LogError($"PolyCubeMapBaker self-test: vertex {v} mask=0x{mask:X2} " +
                               $"classified as {type} (expected Type3).");
                return;
            }
            // Match BakeTexture: the classifying symmetry itself (not its
            // inverse) is what the shader's e.R step needs -- see the comment
            // in BakeTexture for the rationale.
            if (!SymToRotBits.TryGetValue(sym, out var rotBits))
            {
                Debug.LogError($"PolyCubeMapBaker self-test: vertex {v} symmetry not reachable.");
                return;
            }
            try { EncodeByte(type, rotBits); }
            catch (InvalidOperationException ex)
            {
                Debug.LogError($"PolyCubeMapBaker self-test: cannot encode vertex {v}: {ex.Message}");
                return;
            }
            t3Count++;
        }
        if (t3Count != 8)
        {
            Debug.LogError($"PolyCubeMapBaker self-test: expected 8 type-3 cells, got {t3Count}.");
        }
    }
}
