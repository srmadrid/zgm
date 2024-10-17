const std = @import("std");
const zgm = @import("zgm.zig");

pub fn main() void {
    //std.debug.print("Size of @Vector(4, f64): {}\n", .{@sizeOf(@Vector(4, f32))});
    //std.debug.print("Size of zgm.Vector4(f64): {}\n\n", .{@sizeOf(zgm.Vector4(f32))});
    //std.debug.print("Size of [4]@Vector(4, f64): {}\n", .{@sizeOf([4]@Vector(4, f32))});
    //std.debug.print("Size of zgm.Matrix4x4(f64): {}\n", .{@sizeOf(zgm.Matrix4x4(f32))});

    if (false) {
        const n = 100_000_000;
        const A = zgm.Matrix4x4(f32).init(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0);
        const B = zgm.Matrix4x4(f32).init(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0);

        var C: zgm.Matrix4x4(f32) = undefined;

        std.debug.print("{} executions\n", .{n});

        var start_time = std.time.nanoTimestamp();

        for (0..n) |_| {
            C = A.mulMatrix4x4(&B);
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
            C = A.mullMatrix4x4(&B);
        }

        end_time = std.time.nanoTimestamp();
        duration_ns = end_time - start_time;

        // Convert nanoseconds to seconds as a floating-point number.
        // duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));
        duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0);

        // Print the duration in seconds with high precision (e.g., 9 decimal places).
        std.debug.print("Dot mult: {d:.9} seconds\n\n", .{duration_s});

        std.debug.print("1 execution\n", .{});

        start_time = std.time.nanoTimestamp();

        C = A.mulMatrix4x4(&B);

        end_time = std.time.nanoTimestamp();
        duration_ns = end_time - start_time;

        // Convert nanoseconds to seconds as a floating-point number.
        // var duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));
        duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0);

        // Print the duration in seconds with high precision (e.g., 9 decimal places).
        std.debug.print("Unwrapped mult: {d:.9} seconds\n", .{duration_s});

        start_time = std.time.nanoTimestamp();

        C = A.mullMatrix4x4(&B);

        end_time = std.time.nanoTimestamp();
        duration_ns = end_time - start_time;

        // Convert nanoseconds to seconds as a floating-point number.
        // duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));
        duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0);

        // Print the duration in seconds with high precision (e.g., 9 decimal places).
        std.debug.print("Dot mult: {d:.9} seconds\n", .{duration_s});
    }
}
