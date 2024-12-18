const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const module = b.addModule("zgm", .{
        .root_source_file = b.path("src/zgm.zig"),
    });

    // Library
    const lib = b.addStaticLibrary(.{
        .name = "zgm",
        .root_source_file = b.path("src/zgm.zig"),
        .target = target,
        .optimize = optimize,
    });
    b.installArtifact(lib);

    // Executable (for testing)
    const exe = b.addExecutable(.{
        .name = "zgm",
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });
    b.installArtifact(exe);
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    const run_step = b.step("run", "Run the executable");
    run_step.dependOn(&run_cmd.step);

    // Tests
    const lib_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/zgm.zig"),
        .target = target,
        .optimize = optimize,
    });
    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    // Examples
    const example_step = b.step("examples", "Build examples");
    for ([_][]const u8{
        "ndarray_add",
    }) |example_name| {
        const example = b.addExecutable(.{
            .name = example_name,
            .root_source_file = b.path(b.fmt("examples/{s}.zig", .{example_name})),
            .target = target,
            .optimize = optimize,
        });
        const install_example = b.addInstallArtifact(example, .{});
        example.root_module.addImport("zgm", module);
        example_step.dependOn(&example.step);
        example_step.dependOn(&install_example.step);
    }

    // Steps
    const check_step = b.step("check", "Check if the code compiles; this is for ZLS.");
    check_step.dependOn(&exe.step);
}
