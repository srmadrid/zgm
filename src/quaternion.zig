const std = @import("std");
const zgm = @import("zgm.zig");
const Vector3 = zgm.Vector3;
const Vector4 = zgm.Vector4;

/// A quaternion.
pub fn Quaternion(comptime T: type) type {
    return struct {
        v: @Vector(4, T),

        /// Initialize a quaternion with the given components.
        ///
        /// **Parameters**:
        /// - `xs`: The x component.
        /// - `ys`: The y component.
        /// - `zs`: The z component.
        /// - `ws`: The w component.
        ///
        /// **Returns**:
        /// - A new quaternion with the given components.
        pub inline fn init(xs: T, ys: T, zs: T, ws: T) Quaternion(T) {
            return Quaternion(T){
                .v = @Vector(4, T){ xs, ys, zs, ws },
            };
        }

        /// Initialize a quaternion with the given vector.
        ///
        /// **Parameters**:
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - A new quaternion with the given quaternion.
        pub inline fn fromVector(v: *const Vector4(T)) Quaternion(T) {
            return Quaternion(T){
                .v = @Vector(4, T){ v.v[0], v.v[1], v.v[2], v.v[3] },
            };
        }

        /// Initialize a left-handed quaternion from the given axis and angle.
        ///
        /// **Parameters**:
        /// - `axis`: The axis. Must be normalized.
        /// - `angle`: The angle.
        ///
        /// **Returns**:
        /// - A new quaternion from the given axis and angle.
        pub inline fn fromAxisAngleLH(axis: *const Vector3(T), angle: T) Quaternion(T) {
            const halfAngle = angle / 2;
            const sinHalfAngle = @sin(-halfAngle);
            const cosHalfAngle = @cos(halfAngle);

            return Quaternion(T){
                .v = @Vector(4, T){
                    axis.v[0] * sinHalfAngle,
                    axis.v[1] * sinHalfAngle,
                    axis.v[2] * sinHalfAngle,
                    cosHalfAngle,
                },
            };
        }

        /// Initialize a right-handed quaternion from the given axis and angle.
        ///
        /// **Parameters**:
        /// - `axis`: The axis. Must be normalized.
        /// - `angle`: The angle.
        ///
        /// **Returns**:
        /// - A new quaternion from the given axis and angle.
        pub inline fn fromAxisAngleRH(axis: *const Vector3(T), angle: T) Quaternion(T) {
            const halfAngle = angle / 2;
            const sinHalfAngle = @sin(halfAngle);
            const cosHalfAngle = @cos(halfAngle);

            return Quaternion(T){
                .v = @Vector(4, T){
                    axis.v[0] * sinHalfAngle,
                    axis.v[1] * sinHalfAngle,
                    axis.v[2] * sinHalfAngle,
                    cosHalfAngle,
                },
            };
        }

        /// Initialize a quaternion from the given vectors.
        ///
        /// **Parameters**:
        /// - `v`: The first vector. Must be normalized.
        /// - `u`: The second vector. Must be normalized.
        ///
        /// **Returns**:
        /// - A new quaternion from the given vectors.
        pub inline fn fromPair(v: *const Vector3(T), u: *const Vector3(T)) Quaternion(T) {
            const dot = v.dot(u);
            const cross = v.cross(u);

            return Quaternion(T){
                .v = @Vector(4, T){
                    cross.v[0],
                    cross.v[1],
                    cross.v[2],
                    1 + dot,
                },
            };
        }

        /// Initialize a left-handed quaternion from the given Euler angles.
        ///
        /// **Parameters**:
        /// - `pitch`: The pitch angle.
        /// - `yaw`: The yaw angle.
        /// - `roll`: The roll angle.
        ///
        /// **Returns**:
        /// - A new quaternion from the given Euler angles.
        pub inline fn fromEulerLH(pitch: T, yaw: T, roll: T) Quaternion(T) {
            const halfPitch = pitch / 2;
            const sp = @sin(halfPitch);
            const cp = @cos(halfPitch);

            const halfYaw = yaw / 2;
            const sy = @sin(halfYaw);
            const cy = @cos(halfYaw);

            const halfRoll = roll / 2;
            const sr = @sin(halfRoll);
            const cr = @cos(halfRoll);

            return Quaternion(T){
                .v = @Vector(4, T){
                    -(cy * cp * sr - sy * sp * cr),
                    -(cy * sp * cr + sy * cp * sr),
                    -(sy * cp * cr - cy * sp * sr),
                    cy * cp * cr + sy * sp * sr,
                },
            };
        }

        /// Initialize a right-handed quaternion from the given Euler angles.
        ///
        /// **Parameters**:
        /// - `pitch`: The pitch angle.
        /// - `yaw`: The yaw angle.
        /// - `roll`: The roll angle.
        ///
        /// **Returns**:
        /// - A new quaternion from the given Euler angles.
        pub inline fn fromEulerRH(pitch: T, yaw: T, roll: T) Quaternion(T) {
            const halfPitch = pitch / 2;
            const sp = @sin(halfPitch);
            const cp = @cos(halfPitch);

            const halfYaw = yaw / 2;
            const sy = @sin(halfYaw);
            const cy = @cos(halfYaw);

            const halfRoll = roll / 2;
            const sr = @sin(halfRoll);
            const cr = @cos(halfRoll);

            return Quaternion(T){
                .v = @Vector(4, T){
                    cy * cp * sr - sy * sp * cr,
                    cy * sp * cr + sy * cp * sr,
                    sy * cp * cr - cy * sp * sr,
                    cy * cp * cr + sy * sp * sr,
                },
            };
        }

        pub inline fn x(self: *const Quaternion(T)) T {
            return self.v[0];
        }

        pub inline fn y(self: *const Quaternion(T)) T {
            return self.v[1];
        }

        pub inline fn z(self: *const Quaternion(T)) T {
            return self.v[2];
        }

        pub inline fn w(self: *const Quaternion(T)) T {
            return self.v[3];
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - A new quaternion with the element-wise sum of the input
        ///   quaternions.
        pub inline fn add(q: *const Quaternion(T), r: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = q.v + r.v,
            };
        }

        /// Adds a scalar to all components of a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new quaternion with the scalar added to all components.
        pub inline fn addScalar(q: *const Quaternion(T), s: T) Quaternion(T) {
            return Quaternion(T){
                .v = q.v + @Vector(4, T){ s, s, s, s },
            };
        }

        /// Element-wise subtraction.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - A new quaternion with the element-wise difference of the input
        ///   quaternions.
        pub inline fn sub(q: *const Quaternion(T), r: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = q.v - r.v,
            };
        }

        /// Subtracts a scalar from all components of a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new quaternion with the scalar subtracted from all components.
        pub inline fn subScalar(q: *const Quaternion(T), s: T) Vector4(T) {
            return Quaternion(T){
                .v = q.v - @Vector(4, T){ s, s, s, s },
            };
        }

        /// Negates all components of a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        ///
        /// **Returns**:
        /// - A new quaternion with the negated components.
        pub inline fn negate(q: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = -q.v,
            };
        }

        /// Element-wise multiplication.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - A new quaternion with the element-wise product of the input
        ///   quaternions.
        pub inline fn mul(q: *const Quaternion(T), r: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = q.v * r.v,
            };
        }

        /// Multiplies all components of a quaternion by a scalar.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new quaternion with all components multiplied by the scalar.
        pub inline fn mulScalar(q: *const Quaternion(T), s: T) Quaternion(T) {
            return Quaternion(T){
                .v = q.v * @Vector(4, T){ s, s, s, s },
            };
        }
    };
}
