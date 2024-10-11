const std = @import("std");
const zgm = @import("zgm.zig");

pub fn main() void {
    const n = 1_000_000_000;
    const v = zgm.Vector2(f64).init(10, 20);
    const A = zgm.Matrix2x2(f64).init(1, 2, 3, 4);

    var u: zgm.Vector2(f64) = undefined;

    var start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        u = v.mulMatrix2x2(&A);
    }

    var end_time = std.time.nanoTimestamp();
    var duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    // var duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));
    var duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0);

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Unwrapped mult: {d:.9} seconds\n", .{duration_s});

    start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        u = v.mullMatrix2x2(&A);
    }

    end_time = std.time.nanoTimestamp();
    duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    // duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));
    duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0);

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Dot mult: {d:.9} seconds\n", .{duration_s});
}
