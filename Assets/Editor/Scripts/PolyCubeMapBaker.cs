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
