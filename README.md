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
