const std = @import("std");

pub const Vector2 = @import("vector.zig").Vector2;
pub const Vector3 = @import("vector.zig").Vector3;
pub const Vector4 = @import("vector.zig").Vector4;

pub const Matrix2x2 = @import("matrix.zig").Matrix2x2;
pub const Matrix2x3 = @import("matrix.zig").Matrix2x3;
pub const Matrix2x4 = @import("matrix.zig").Matrix2x4;
pub const Matrix3x2 = @import("matrix.zig").Matrix3x2;
pub const Matrix3x3 = @import("matrix.zig").Matrix3x3;
pub const Matrix3x4 = @import("matrix.zig").Matrix3x4;
pub const Matrix4x2 = @import("matrix.zig").Matrix4x2;
pub const Matrix4x3 = @import("matrix.zig").Matrix4x3;
pub const Matrix4x4 = @import("matrix.zig").Matrix4x4;

pub const Quaternion = @import("quaternion.zig").Quaternion;

pub const RGB = @import("color.zig").RGB;
pub const RGBA = @import("color.zig").RGBA;
pub const HSL = @import("color.zig").HSL;
pub const HSLA = @import("color.zig").HSLA;
pub const HSV = @import("color.zig").HSV;
pub const HSVA = @import("color.zig").HSVA;
pub const CMYK = @import("color.zig").CMYK;

test {
    std.testing.refAllDeclsRecursive(@This());
}
