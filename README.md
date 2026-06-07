Unity implementation of PolyCubeMaps (original code by Marco Tarini):<br>
https://vcg.isti.cnr.it/polycubemaps/

![image](https://github.com/user-attachments/assets/7084ffa8-8936-48a3-a27c-22bcc4439e57)

## Generating PolyCubeMaps

The repo also includes a sample baker that demonstrates how to *generate*
PolyCubeMap data (3D UVs + LUT + texture) from a mesh, following
Tarini et al. 2004. The implementation lives in
`Assets/Editor/Scripts/PolyCubeMapBaker.cs` and exposes two Editor menu
items:

- **PolyCubeMap → Bake → Sphere Sample (single-cube polycube)** — bakes
  a procedurally-generated icosphere onto a single-cube polycube (the
  "apple" example from paper §3, Figure 6: 8 dual cells, all type 3).
- **PolyCubeMap → Bake → L-shape Sample** — bakes an L-shaped polycube
  (3 cubes) with the polycube's own surface as the mesh; exercises
  cell types 3, 4a, and 4b.

Each menu item:

1. Computes per-vertex 3D texture coordinates by warping the mesh AABB
   into the polycube AABB and projecting each vertex to the closest
   point on the polycube surface (paper §6.1, simplified — the iterative
   extended-MIPS optimization is left as a documented stub).
2. Classifies every polycube vertex into a cell type (3 / 4a / 4b / 5)
   and a rotation in the 5-bit op space the shader decodes (paper §5.1).
3. Bakes a `1024×1024` RGBA32 texture containing the T3LUT and one
   square-packed patch per non-empty cell, filled with a per-cell
   HSV-hashed checkerboard so seamlessness and patch placement are
   visible at a glance.
4. Writes the mesh and texture to `Assets/Samples/Generated/`, creates a
   `GameObject` in the active scene with the existing `PolyCubeMap`
   component wired up, and logs a bake summary.

Cell types 6a/6b and the iterative parameterization optimization are
intentionally unimplemented — they match the limits of the runtime
shader and are flagged with `TODO` comments and warnings in the baker.

### From Mesh… (arbitrary mesh + optional texture)

A third menu item, **PolyCubeMap → Bake → From Mesh…**, opens an editor
window that lets you bake a PolyCubeMap from any Unity `Mesh` plus an
optional source `Texture2D`. The pipeline:

1. **Voxelize** the mesh on a `resolution × resolution × resolution`
   grid (longest mesh-AABB axis maps to `resolution`). Shell voxels are
   marked by barycentric-stepping every triangle at sub-voxel density;
   if "Fill interior" is on, a 3D BFS flood from outside the bounding
   box marks every voxel it cannot reach as interior. The polycube is
   `shell ∪ interior`. Resolution is hard-capped at 14 because the
   shader packs the LUT as `cellX + 16 * cellZ`.
2. **3D UVs** are computed exactly as for the sample meshes (warp +
   closest-point on the polycube surface).
3. **Bake the atlas**. If a source texture is provided (and is
   Read/Write enabled), each atlas pixel inside a patch is back-projected
   into mesh space, snapped to the nearest mesh-surface point, and
   sampled from the source via the interpolated UV0. Otherwise the
   existing HSV checkerboard fill is used. Meshes with more than 20k
   triangles fall back to a vertex-only nearest search for speed.
4. The resulting mesh + PNG land in `Assets/Samples/Generated/`, and a
   wired-up `GameObject` is spawned and selected in the scene.

This voxelization is a deliberately simple demo — for production-quality
polycube construction, see the Fu/Bai/Liu 2016 algorithm in
[`Material/paper.pdf`](Material/paper.pdf) and the reference C++ in
`Material/*.cpp`. The squarelet-to-facelet projection used by the
texture sampler is also a single uniform scheme rather than the exact
per-case packing the shader assumes, so patches will look approximately
right on the four supported types but not pixel-perfect.
