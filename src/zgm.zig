const std = @import("std");

pub const Vector2 = @import("vector.zig").Vector2;
pub const Vector3 = @import("vector.zig").Vector3;
pub const Vector4 = @import("vector.zig").Vector4;

test {
    std.testing.refAllDeclsRecursive(@This());
}
