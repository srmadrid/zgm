const std = @import("std");
const zgm = @import("zgm.zig");
const Vector3 = zgm.Vector3;
const Vector4 = zgm.Vector4;
const Matrix4x4 = zgm.Matrix4x4;

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

        /// Initialize a quaternion with the given 3D vector.
        ///
        /// **Parameters**:
        /// - `v`: The vector.
        /// - `ws`: The w component.
        ///
        /// **Returns**:
        /// - A new quaternion with the given vector.
        pub inline fn fromVector3(v: *const Vector3(T), ws: T) Quaternion(T) {
            return Quaternion(T){
                .v = @Vector(4, T){ v.v[0], v.v[1], v.v[2], ws },
            };
        }

        /// Initialize a quaternion with the given 4D vector.
        ///
        /// **Parameters**:
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - A new quaternion with the given quaternion.
        pub inline fn fromVector4(v: *const Vector4(T)) Quaternion(T) {
            return Quaternion(T){
                .v = @Vector(4, T){ v.v[0], v.v[1], v.v[2], v.v[3] },
            };
        }

        /// Initialize a quaternion from the given axis and angle.
        ///
        /// **Parameters**:
        /// - `axis`: The axis. Must be normalized.
        /// - `angle`: The angle.
        ///
        /// **Returns**:
        /// - A new quaternion from the given axis and angle.
        pub inline fn fromAxisAngle(axis: *const Vector3(T), angle: T) Quaternion(T) {
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

        /// Initialize a quaternion from rotation needed to rotate one vector to
        /// another.
        ///
        /// **Parameters**:
        /// - `from`: The first vector. Must be normalized.
        /// - `to`: The second vector. Must be normalized.
        ///
        /// **Returns**:
        /// - A new quaternion from the given vectors.
        pub inline fn fromPair(from: *const Vector3(T), to: *const Vector3(T)) Quaternion(T) {
            const dott = from.dot(to);
            const cross = from.cross(to);

            const q = Quaternion(T){
                .v = @Vector(4, T){
                    cross.v[0],
                    cross.v[1],
                    cross.v[2],
                    1 + dott,
                },
            };

            return q.normalize();
        }

        /// Initialize a quaternion from the given Euler angles.
        ///
        /// **Parameters**:
        /// - `pitch`: The pitch angle.
        /// - `yaw`: The yaw angle.
        /// - `roll`: The roll angle.
        ///
        /// **Returns**:
        /// - A new quaternion from the given Euler angles.
        pub inline fn fromEuler(pitch: T, yaw: T, roll: T) Quaternion(T) {
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

        /// Initialize a quaternion from the given rotation matrix.
        ///
        /// **Parameters**:
        /// - `A`: The rotation matrix.
        ///
        /// **Returns**:
        /// - A new quaternion from the given rotation matrix.
        pub inline fn fromMatrix4x4(A: *const Matrix4x4(T)) Quaternion(T) {
            const trace = A.m[0].v[0] + A.m[1].v[1] + A.m[2].v[2];
            const tracePlusOne = trace + 1;

            if (tracePlusOne > 0) {
                const s = 0.5 / @sqrt(tracePlusOne);
                return Quaternion(T){
                    .v = @Vector(4, T){
                        (A.m[1].v[2] - A.m[2].v[1]) * s,
                        (A.m[2].v[0] - A.m[0].v[2]) * s,
                        (A.m[0].v[1] - A.m[1].v[0]) * s,
                        0.25 / s,
                    },
                };
            } else if (A.m[0].v[0] > A.m[1].v[1] and A.m[0].v[0] > A.m[2].v[2]) {
                const s = 2 * @sqrt(1 + A.m[0].v[0] - A.m[1].v[1] - A.m[2].v[2]);
                return Quaternion(T){
                    .v = @Vector(4, T){
                        0.25 * s,
                        (A.m[1].v[0] + A.m[0].v[1]) / s,
                        (A.m[2].v[0] + A.m[0].v[2]) / s,
                        (A.m[1].v[2] - A.m[2].v[1]) / s,
                    },
                };
            } else if (A.m[1].v[1] > A.m[2].v[2]) {
                const s = 2 * @sqrt(1 + A.m[1].v[1] - A.m[0].v[0] - A.m[2].v[2]);
                return Quaternion(T){
                    .v = @Vector(4, T){
                        (A.m[1].v[0] + A.m[0].v[1]) / s,
                        0.25 * s,
                        (A.m[2].v[1] + A.m[1].v[2]) / s,
                        (A.m[2].v[0] - A.m[0].v[2]) / s,
                    },
                };
            } else {
                const s = 2 * @sqrt(1 + A.m[2].v[2] - A.m[0].v[0] - A.m[1].v[1]);
                return Quaternion(T){
                    .v = @Vector(4, T){
                        (A.m[2].v[0] + A.m[0].v[2]) / s,
                        (A.m[2].v[1] + A.m[1].v[2]) / s,
                        0.25 * s,
                        (A.m[0].v[1] - A.m[1].v[0]) / s,
                    },
                };
            }
        }

        /// Initialize a 4x4 rotation matrix from the quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion. Must be normalized.
        ///
        /// **Returns**:
        /// - A new 4x4 rotation matrix from the quaternion.
        pub inline fn toMatrix4x4(q: *const Quaternion(T)) Matrix4x4(T) {
            const xx = q.v[0] * q.v[0];
            const xy = q.v[0] * q.v[1];
            const xz = q.v[0] * q.v[2];
            const xw = q.v[0] * q.v[3];

            const yy = q.v[1] * q.v[1];
            const yz = q.v[1] * q.v[2];
            const yw = q.v[1] * q.v[3];

            const zz = q.v[2] * q.v[2];
            const zw = q.v[2] * q.v[3];

            return Matrix4x4(T){
                .m = [4]Vector4(T){
                    Vector4(T).init(1 - 2 * (yy + zz), 2 * (xy + zw), 2 * (xz - yw), 0),
                    Vector4(T).init(2 * (xy - zw), 1 - 2 * (xx + zz), 2 * (yz + xw), 0),
                    Vector4(T).init(2 * (xz + yw), 2 * (yz - xw), 1 - 2 * (xx + yy), 0),
                    Vector4(T).init(0, 0, 0, 1),
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
        pub inline fn subScalar(q: *const Quaternion(T), s: T) Quaternion(T) {
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

        /// Quaternion-Vector3 multiplication.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The rotated vector.
        pub inline fn mulVector3(q: *const Quaternion(T), v: *const Vector3(T)) Vector3(T) {
            const qv = Quaternion(T).fromVector3(v, 0);
            const qvq = q.mulQuaternion(&qv).mulQuaternion(&q.conjugate());
            return Vector3(T){
                .v = @Vector(3, T){ qvq.v[0], qvq.v[1], qvq.v[2] },
            };
        }

        /// Quaternion-Vector4 multiplication.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The rotated vector.
        pub inline fn mulVector4(q: *const Quaternion(T), v: *const Vector4(T)) Vector4(T) {
            const qv = Quaternion(T).fromVector4(v);
            const qvq = q.mulQuaternion(&qv).mulQuaternion(&q.conjugate());
            return Vector4(T){
                .v = @Vector(4, T){ qvq.v[0], qvq.v[1], qvq.v[2], v.v[3] },
            };
        }

        /// Quaternion-Quaternion multiplication.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - The product of the input quaternions.
        pub inline fn mulQuaternion(q: *const Quaternion(T), r: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = @Vector(4, T){
                    q.v[3] * r.v[0] + q.v[0] * r.v[3] + q.v[1] * r.v[2] - q.v[2] * r.v[1],
                    q.v[3] * r.v[1] - q.v[0] * r.v[2] + q.v[1] * r.v[3] + q.v[2] * r.v[0],
                    q.v[3] * r.v[2] + q.v[0] * r.v[1] - q.v[1] * r.v[0] + q.v[2] * r.v[3],
                    q.v[3] * r.v[3] - q.v[0] * r.v[0] - q.v[1] * r.v[1] - q.v[2] * r.v[2],
                },
            };
        }

        /// Element-wise division.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - A new quaternion with the element-wise division of the input
        ///   quaternions.
        pub inline fn div(q: *const Quaternion(T), r: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = q.v / r.v,
            };
        }

        /// Divides all components of a quaternion by a scalar.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new quaternion with all components divided by the scalar.
        pub inline fn divScalar(q: *const Quaternion(T), s: T) Quaternion(T) {
            return Quaternion(T){
                .v = q.v / @Vector(4, T){ s, s, s, s },
            };
        }

        /// Conjugates a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        ///
        /// **Returns**:
        /// - The conjugate of the input quaternion.
        pub inline fn conjugate(q: *const Quaternion(T)) Quaternion(T) {
            return Quaternion(T){
                .v = @Vector(4, T){ -q.v[0], -q.v[1], -q.v[2], q.v[3] },
            };
        }

        /// Computes the dot product of two quaternions.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - The dot product of the input quaternions.
        pub inline fn dot(q: *const Quaternion(T), r: *const Quaternion(T)) T {
            return @reduce(.Add, q.v * r.v);
        }

        /// Computes the modulus of a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        ///
        /// **Returns**:
        /// - The modulus of the input quaternion.
        pub inline fn mod(q: *const Quaternion(T)) T {
            return @sqrt(@reduce(.Add, q.v * q.v));
        }

        /// Computes the squared modulus of a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        ///
        /// **Returns**:
        /// - The squared modulus of the input quaternion.
        pub inline fn mod2(q: *const Quaternion(T)) T {
            return @reduce(.Add, q.v * q.v);
        }

        /// Inverts a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        ///
        /// **Returns**:
        /// - The inverse of the input quaternion.
        pub inline fn inverse(q: *const Quaternion(T)) Quaternion(T) {
            return q.conjugate().divScalar(q.mod2());
        }

        /// Normalizes a quaternion.
        ///
        /// **Parameters**:
        /// - `q`: The quaternion.
        ///
        /// **Returns**:
        /// - The normalized quaternion.
        pub inline fn normalize(q: *const Quaternion(T)) Quaternion(T) {
            return q.divScalar(q.mod());
        }

        /// Linearly interpolates between two quaternions.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        /// - `t`: The interpolation factor.
        ///
        /// **Returns**:
        /// - The interpolated quaternion.
        pub inline fn lerp(q: *const Quaternion(T), r: *const Quaternion(T), t: T) Quaternion(T) {
            return q.mulScalar(1 - t).add(&r.mulScalar(t));
        }

        /// Normalized linear interpolation between two quaternions.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        /// - `t`: The interpolation factor.
        ///
        /// **Returns**:
        /// - The normalized interpolated quaternion.
        pub inline fn nlerp(q: *const Quaternion(T), r: *const Quaternion(T), t: T) Quaternion(T) {
            return q.lerp(r, t).normalize();
        }

        /// Spherical linear interpolation between two quaternions.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        /// - `t`: The interpolation factor.
        ///
        /// **Returns**:
        /// - The interpolated quaternion.
        pub inline fn slerp(q: *const Quaternion(T), r: *const Quaternion(T), t: T) Quaternion(T) {
            const dott = q.dot(r);
            const threshold = 0.9995;

            if (dott > threshold) {
                return q.nlerp(r, t);
            }

            const angle = std.math.acos(dott);
            const invSinAngle = 1 / @sin(angle);
            const qFactor = @sin((1 - t) * angle) * invSinAngle;
            const rFactor = @sin(t * angle) * invSinAngle;

            return q.mulScalar(qFactor).add(&r.mulScalar(rFactor));
        }

        /// Compares two quaternions and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        ///
        /// **Returns**:
        /// - `true` if the quaternions are equal, `false` otherwise.
        pub inline fn equal(q: *const Quaternion(T), r: *const Quaternion(T)) bool {
            return q.v[0] == r.v[0] and q.v[1] == r.v[1] and q.v[2] == r.v[2] and q.v[3] == r.v[3];
        }

        /// Compares two quaternions and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `q`: The first quaternion.
        /// - `r`: The second quaternion.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the quaternions are approximately equal, `false`
        ///   otherwise.
        pub inline fn approxEqual(q: *const Quaternion(T), r: *const Quaternion(T), tolerance: T) bool {
            return @abs(q.v[0] - r.v[0]) <= tolerance and @abs(q.v[1] - r.v[1]) <= tolerance and @abs(q.v[2] - r.v[2]) <= tolerance and @abs(q.v[3] - r.v[3]) <= tolerance;
        }
    };
}

test "Quaternion.add" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.add(&q2);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(6, 8, 10, 12), 0.0001));
}

test "Quaternion.addScalar" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.addScalar(5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(6, 7, 8, 9), 0.0001));
}

test "Quaternion.sub" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.sub(&q2);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(-4, -4, -4, -4), 0.0001));
}

test "Quaternion.subScalar" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.subScalar(5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(-4, -3, -2, -1), 0.0001));
}

test "Quaternion.negate" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.negate();

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(-1, -2, -3, -4), 0.0001));
}

test "Quaternion.mul" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.mul(&q2);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(5, 12, 21, 32), 0.0001));
}

test "Quaternion.mulScalar" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.mulScalar(5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(5, 10, 15, 20), 0.0001));
}

test "Quaternion.mulVector3" {
    const q = Quaternion(f32).fromAxisAngle(&Vector3(f32).init(1, 0, 0).normalize(), std.math.pi / 2.0);
    const v = Vector3(f32).init(0, 1, 0);
    const r = q.mulVector3(&v);

    try std.testing.expect(r.approxEqual(&Vector3(f32).init(0, 0, 1), 0.0001));
}

test "Quaternion.mulVector4" {
    const q = Quaternion(f32).fromAxisAngle(&Vector3(f32).init(1, 0, 0).normalize(), std.math.pi / 2.0);
    const v = Vector4(f32).init(0, 1, 0, 1);
    const r = q.mulVector4(&v);

    try std.testing.expect(r.approxEqual(&Vector4(f32).init(0, 0, 1, 1), 0.0001));
}

test "Quaternion.mulQuaternion" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.mulQuaternion(&q2);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(24, 48, 48, -6), 0.0001));
}

test "Quaternion.div" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.div(&q2);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(0.2, 0.3333, 0.4285, 0.5), 0.0001));
}

test "Quaternion.divScalar" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.divScalar(5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(0.2, 0.4, 0.6, 0.8), 0.0001));
}

test "Quaternion.conjugate" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.conjugate();

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(-1, -2, -3, 4), 0.0001));
}

test "Quaternion.dot" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.dot(&q2);

    try std.testing.expectApproxEqRel(70, r, 0.001);
}

test "Quaternion.mod" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.mod();

    try std.testing.expectApproxEqRel(5.4772, r, 0.0001);
}

test "Quaternion.mod2" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.mod2();

    try std.testing.expectApproxEqRel(30, r, 0.0001);
}

test "Quaternion.inverse" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.inverse();

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(-0.03333, -0.06666, -0.1, 0.1333), 0.0001));
}

test "Quaternion.normalize" {
    const q = Quaternion(f32).init(1, 2, 3, 4);
    const r = q.normalize();

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(0.1826, 0.3651, 0.5477, 0.7302), 0.0001));
}

test "Quaternion.lerp" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.lerp(&q2, 0.5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(3, 4, 5, 6), 0.0001));
}

test "Quaternion.nlerp" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.nlerp(&q2, 0.5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(0.3234, 0.4313, 0.5391, 0.6469), 0.0001));
}

test "Quaternion.slerp" {
    const q1 = Quaternion(f32).init(1, 2, 3, 4);
    const q2 = Quaternion(f32).init(5, 6, 7, 8);
    const r = q1.slerp(&q2, 0.5);

    try std.testing.expect(r.approxEqual(&Quaternion(f32).init(0.3234, 0.4313, 0.5391, 0.6469), 0.0001));
}

test "Quaternion.fromAxisAngle" {
    const q = Quaternion(f32).fromAxisAngle(&Vector3(f32).init(1, 0, 0), std.math.pi / 2.0);
    const v = Vector3(f32).init(0, 1, 0);
    const r = q.mulVector3(&v);

    try std.testing.expect(r.approxEqual(&Vector3(f32).init(0, 0, 1), 0.0001));
}

test "Quaternion.fromPair" {
    const q = Quaternion(f32).fromPair(&Vector3(f32).init(1, 0, 0).normalize(), &Vector3(f32).init(0, 1, 0).normalize());
    const v = Vector3(f32).init(0, 1, 0);
    const r = q.mulVector3(&v);

    try std.testing.expect(r.approxEqual(&Vector3(f32).init(-1, 0, 0), 0.0001));
}

test "Quaternion.fromEuler" {
    const q = Quaternion(f32).fromEuler(std.math.pi / 2.0, 0, 0);
    const v = Vector3(f32).init(0, 1, 0);
    const r = q.mulVector3(&v);

    try std.testing.expect(r.approxEqual(&Vector3(f32).init(0, 1, 0), 0.0001));
}

test "Quaternion.fromMatrix4x4" {
    const q = Quaternion(f32).fromMatrix4x4(&Matrix4x4(f32).rotate(std.math.pi / 2.0, &Vector3(f32).init(1, 0, 0)));
    const v = Vector3(f32).init(0, 1, 0);
    const r = q.mulVector3(&v);

    try std.testing.expect(r.approxEqual(&Vector3(f32).init(0, 0, 1), 0.0001));
}

test "Quaternion.toMatrix4x4" {
    const q = Quaternion(f32).fromAxisAngle(&Vector3(f32).init(1, 0, 0), std.math.pi / 2.0);
    const r = q.toMatrix4x4();

    try std.testing.expect(r.approxEqual(&Matrix4x4(f32).rotate(std.math.pi / 2.0, &Vector3(f32).init(1, 0, 0)), 0.0001));
}
