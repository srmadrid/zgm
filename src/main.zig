const std = @import("std");

const zgm = @import("zgm.zig");

pub fn main() void {
    const v = zgm.Vector2(f64).init(10, 20);
    const w = zgm.Vector2(f64).init(20, 40);

    const x = zgm.Vector2(f64).add(&v, &w);

    std.debug.print("{}\n", .{x.x()});
}
