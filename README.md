# ZGM: Zig Graphics Mathematics

A Zig library for mathematics for graphics programming.

## Current Features

- Vectors
- Matrices
  - Column-major

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
