# ZGM: Zig Graphics Mathematics

A Zig library for mathematics for graphics programming.

## Current Features

- Vectors
- Matrices
  - Column-major
  - Transformations are offered only for homogeneous coordinates, i.e., 3x3 matrices work as 2D transformations and 4x4 matrices work as 3D transformations, and no 2x2 or 3x3 matrices are provided for 2D or 3D transformations, respectively.
- Quaternions
- Colors

All transformations assume a right-handed coordinate system.

## Installation

To use this library in your project, run

```bash
zig fetch --save git+https://github.com/srmadrid/zgm
```

and add it to your `build.zig` file:

```zig
const zgm = b.dependency("zgm", .{});
exe.root_module.addImport("zgm", zgm.module("zgm"));
```
