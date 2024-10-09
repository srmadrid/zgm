const std = @import("std");

/// A two-dimensional vector.
pub fn Vector2(comptime T: type) type {
    return struct {
        v: @Vector(2, T),

        /// Initialize a vector with the given components.
        pub inline fn init(xs: T, ys: T) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ xs, ys },
            };
        }

        pub inline fn x(self: *const Vector2(T)) T {
            return self.v[0];
        }

        pub inline fn y(self: *const Vector2(T)) T {
            return self.v[1];
        }

        /// Splat a scalar into a vector.
        pub inline fn splat(s: T) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ s, s },
            };
        }

        /// Element-wise addition.
        pub inline fn add(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = v.v + w.v,
            };
        }

        /// Scalar addition.
        pub inline fn addScalar(v: *const Vector2(T), s: T) Vector2(T) {
            return Vector2(T){
                .v = v.v + @Vector(2, T){ s, s },
            };
        }

        /// Element-wise subtraction.
        pub inline fn sub(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = v.v - w.v,
            };
        }

        /// Scalar subtraction.
        pub inline fn subScalar(v: *const Vector2(T), s: T) Vector2(T) {
            return Vector2(T){
                .v = v.v - @Vector(2, T){ s, s },
            };
        }

        /// Negates a vector.
        pub inline fn negate(v: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ -v.v[0], -v.v[1] },
            };
        }

        /// Element-wise multiplication.
        pub inline fn mul(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = v.v * w.v,
            };
        }

        /// Scalar multiplication.
        pub inline fn mulScalar(v: *const Vector2(T), s: T) Vector2(T) {
            return Vector2(T){
                .v = v.v * @Vector(2, T){ s, s },
            };
        }

        /// Element-wise division.
        pub inline fn div(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = v.v / w.v,
            };
        }

        /// Scalar division.
        pub inline fn divScalar(v: *const Vector2(T), s: T) Vector2(T) {
            return Vector2(T){
                .v = v.v / @Vector(2, T){ s, s },
            };
        }

        /// Computes the modulus of a vector.
        pub inline fn mod(v: *const Vector2(T)) T {
            return @sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1]);
        }

        /// Computes the squared modulus of a vector.
        pub inline fn mod2(v: *const Vector2(T)) T {
            return v.v[0] * v.v[0] + v.v[1] * v.v[1];
        }

        /// Normalizes a vector.
        pub inline fn normalize(v: *const Vector2(T)) Vector2(T) {
            return v.div(&splat(v.mod()));
        }

        /// Computes the dot product between two vectors.
        pub inline fn dot(v: *const Vector2(T), w: *const Vector2(T)) T {
            return v.v[0] * w.v[0] + v.v[1] * w.v[1];
        }

        /// Computes the distance between two vectors.
        pub inline fn distance(v: *const Vector2(T), w: *const Vector2(T)) T {
            return v.sub(w).mod();
        }

        /// Computes the squared distance between two vectors.
        pub inline fn distance2(v: *const Vector2(T), w: *const Vector2(T)) T {
            return v.sub(w).mod2();
        }

        /// Linearly interpolates two vectors.
        pub inline fn lerp(v: *const Vector2(T), w: *const Vector2(T), t: T) Vector2(T) {
            return v.mulScalar(1 - t).add(&w.mulScalar(t));
        }

        /// Clamps element-wise a vector between a maximum and a minimum.
        pub inline fn clamp(v: *const Vector2(T), mn: T, mx: T) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ @min(@max(v.v[0], mn), mx), @min(@max(v.v[1], mn), mx) },
            };
        }

        /// Computes the angle between two vectors.
        pub inline fn angle(v: *const Vector2(T), w: *const Vector2(T)) T {
            return std.math.acos((v.dot(w)) / (v.mod() * w.mod()));
        }

        /// Computes the projection of the first vector onto the second.
        pub inline fn project(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return w.mul(&splat(v.dot(w) / w.mod2()));
        }

        /// Reflects a vector against a normal.
        pub inline fn reflect(v: *const Vector2(T), normal: *const Vector2(T)) Vector2(T) {
            return v.sub(&normal.mul(&splat(2 * v.dot(normal) / normal.mod2())));
        }

        /// Refracts a vector against a normal.
        pub inline fn refract(v: *const Vector2(T), normal: *const Vector2(T), eta: T) Vector2(T) {
            const v_norm = v.normalize();
            const k = 1 - eta * eta * (1 - v_norm.dot(normal) * v_norm.dot(normal));

            if (k < 0) {
                return splat(0);
            }

            return v_norm.mulScalar(eta).sub(&normal.mulScalar(eta * v_norm.dot(normal) + @sqrt(k)));
        }

        /// Returns a vector with the maximum values between two vectors.
        pub inline fn max(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ @max(v.v[0], w.v[0]), @max(v.v[1], w.v[1]) },
            };
        }

        /// Returns a vector with the minimum values between two vectors.
        pub inline fn min(v: *const Vector2(T), w: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ @min(v.v[0], w.v[0]), @min(v.v[1], w.v[1]) },
            };
        }

        /// Returns a vector with the absolute values of the input vector.
        pub inline fn abs(v: *const Vector2(T)) Vector2(T) {
            return Vector2(T){
                .v = @Vector(2, T){ @abs(v.v[0]), @abs(v.v[1]) },
            };
        }

        /// Compares two vectors and returns true if they are equal.
        pub inline fn equal(v: *const Vector2(T), w: *const Vector2(T)) bool {
            return v.v[0] == w.v[0] and v.v[1] == w.v[1];
        }

        /// Compares two vectors and returns true if they are approximately equal.
        pub inline fn approxEqual(v: *const Vector2(T), w: *const Vector2(T), tolerance: T) bool {
            return @abs(v.v[0] - w.v[0]) <= tolerance and @abs(v.v[1] - w.v[1]) <= tolerance;
        }
    };
}

/// A three-dimensional vector.
pub fn Vector3(comptime T: type) type {
    return struct {
        v: @Vector(3, T),

        /// Initialize a vector with the given components.
        pub inline fn init(xs: T, ys: T, zs: T) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ xs, ys, zs },
            };
        }

        pub inline fn x(self: *const Vector3(T)) T {
            return self.v[0];
        }

        pub inline fn y(self: *const Vector3(T)) T {
            return self.v[1];
        }

        pub inline fn z(self: *const Vector3(T)) T {
            return self.v[2];
        }

        /// Splat a scalar into a vector.
        pub inline fn splat(s: T) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ s, s, s },
            };
        }

        /// Element-wise addition.
        pub inline fn add(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = v.v + w.v,
            };
        }

        /// Scalar addition.
        pub inline fn addScalar(v: *const Vector3(T), s: T) Vector3(T) {
            return Vector3(T){
                .v = v.v + @Vector(3, T){ s, s, s },
            };
        }

        /// Element-wise subtraction.
        pub inline fn sub(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = v.v - w.v,
            };
        }

        /// Scalar subtraction.
        pub inline fn subScalar(v: *const Vector3(T), s: T) Vector3(T) {
            return Vector3(T){
                .v = v.v - @Vector(3, T){ s, s, s },
            };
        }

        /// Negates a vector.
        pub inline fn negate(v: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ -v.v[0], -v.v[1], -v.v[2] },
            };
        }

        /// Element-wise multiplication.
        pub inline fn mul(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = v.v * w.v,
            };
        }

        /// Scalar multiplication.
        pub inline fn mulScalar(v: *const Vector3(T), s: T) Vector3(T) {
            return Vector3(T){
                .v = v.v * @Vector(3, T){ s, s, s },
            };
        }

        /// Element-wise division.
        pub inline fn div(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = v.v / w.v,
            };
        }

        /// Scalar division.
        pub inline fn divScalar(v: *const Vector3(T), s: T) Vector3(T) {
            return Vector3(T){
                .v = v.v / @Vector(3, T){ s, s, s },
            };
        }

        /// Computes the modulus of a vector.
        pub inline fn mod(v: *const Vector3(T)) T {
            return @sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]);
        }

        /// Computes the squared modulus of a vector.
        pub inline fn mod2(v: *const Vector3(T)) T {
            return v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2];
        }

        /// Normalizes a vector.
        pub inline fn normalize(v: *const Vector3(T)) Vector3(T) {
            return v.div(&splat(v.mod()));
        }

        /// Computes the dot product between two vectors.
        pub inline fn dot(v: *const Vector3(T), w: *const Vector3(T)) T {
            return v.v[0] * w.v[0] + v.v[1] * w.v[1] + v.v[2] * w.v[2];
        }

        /// Computes the cross product between two vectors.
        pub inline fn cross(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){
                    v.v[1] * w.v[2] - v.v[2] * w.v[1],
                    v.v[2] * w.v[0] - v.v[0] * w.v[2],
                    v.v[0] * w.v[1] - v.v[1] * w.v[0],
                },
            };
        }

        /// Computes the distance between two vectors.
        pub inline fn distance(v: *const Vector3(T), w: *const Vector3(T)) T {
            return v.sub(w).mod();
        }

        /// Computes the squared distance between two vectors.
        pub inline fn distance2(v: *const Vector3(T), w: *const Vector3(T)) T {
            return v.sub(w).mod2();
        }

        /// Linearly interpolates two vectors.
        pub inline fn lerp(v: *const Vector3(T), w: *const Vector3(T), t: T) Vector3(T) {
            return v.mulScalar(1 - t).add(&w.mulScalar(t));
        }

        /// Clamps element-wise a vector between a maximum and a minimum.
        pub inline fn clamp(v: *const Vector3(T), mn: T, mx: T) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ @min(@max(v.v[0], mn), mx), @min(@max(v.v[1], mn), mx), @min(@max(v.v[2], mn), mx) },
            };
        }

        /// Computes the angle between two vectors.
        pub inline fn angle(v: *const Vector3(T), w: *const Vector3(T)) T {
            return std.math.acos((v.dot(w)) / (v.mod() * w.mod()));
        }

        /// Computes the projection of the first vector onto the second.
        pub inline fn project(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return w.mul(&splat(v.dot(w) / w.mod2()));
        }

        /// Reflects a vector against a normal.
        pub inline fn reflect(v: *const Vector3(T), normal: *const Vector3(T)) Vector3(T) {
            return v.sub(&normal.mul(&splat(2 * v.dot(normal) / normal.mod2())));
        }

        /// Refracts a vector against a normal.
        pub inline fn refract(v: *const Vector3(T), normal: *const Vector3(T), eta: T) Vector3(T) {
            const v_norm = v.normalize();
            const k = 1 - eta * eta * (1 - v_norm.dot(normal) * v_norm.dot(normal));

            if (k < 0) {
                return splat(0);
            }

            return v_norm.mulScalar(eta).sub(&normal.mulScalar(eta * v_norm.dot(normal) + @sqrt(k)));
        }

        /// Returns a vector with the maximum values between two vectors.
        pub inline fn max(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ @max(v.v[0], w.v[0]), @max(v.v[1], w.v[1]), @max(v.v[2], w.v[2]) },
            };
        }

        /// Returns a vector with the minimum values between two vectors.
        pub inline fn min(v: *const Vector3(T), w: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ @min(v.v[0], w.v[0]), @min(v.v[1], w.v[1]), @min(v.v[2], w.v[2]) },
            };
        }

        /// Returns a vector with the absolute values of the input vector.
        pub inline fn abs(v: *const Vector3(T)) Vector3(T) {
            return Vector3(T){
                .v = @Vector(3, T){ @abs(v.v[0]), @abs(v.v[1]), @abs(v.v[2]) },
            };
        }

        /// Compares two vectors and returns true if they are equal.
        pub inline fn equal(v: *const Vector3(T), w: *const Vector3(T)) bool {
            return v.v[0] == w.v[0] and v.v[1] == w.v[1] and v.v[2] == w.v[2];
        }

        /// Compares two vectors and returns true if they are approximately equal.
        pub inline fn approxEqual(v: *const Vector3(T), w: *const Vector3(T), tolerance: T) bool {
            return @abs(v.v[0] - w.v[0]) <= tolerance and @abs(v.v[1] - w.v[1]) <= tolerance and @abs(v.v[2] - w.v[2]) <= tolerance;
        }
    };
}

/// A four-dimensional vector.
pub fn Vector4(comptime T: type) type {
    return struct {
        v: @Vector(4, T),

        /// Initialize a vector with the given components.
        pub inline fn init(xs: T, ys: T, zs: T, ws: T) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ xs, ys, zs, ws },
            };
        }

        pub inline fn x(self: *const Vector4(T)) T {
            return self.v[0];
        }

        pub inline fn y(self: *const Vector4(T)) T {
            return self.v[1];
        }

        pub inline fn z(self: *const Vector4(T)) T {
            return self.v[2];
        }

        pub inline fn w(self: *const Vector4(T)) T {
            return self.v[3];
        }

        /// Splat a scalar into a vector.
        pub inline fn splat(s: T) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ s, s, s, s },
            };
        }

        /// Element-wise addition.
        pub inline fn add(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = v.v + u.v,
            };
        }

        /// Scalar addition.
        pub inline fn addScalar(v: *const Vector4(T), s: T) Vector4(T) {
            return Vector4(T){
                .v = v.v + @Vector(4, T){ s, s, s, s },
            };
        }

        /// Element-wise subtraction.
        pub inline fn sub(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = v.v - u.v,
            };
        }

        /// Scalar subtraction.
        pub inline fn subScalar(v: *const Vector4(T), s: T) Vector4(T) {
            return Vector4(T){
                .v = v.v - @Vector(4, T){ s, s, s, s },
            };
        }

        /// Negates a vector.
        pub inline fn negate(v: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ -v.v[0], -v.v[1], -v.v[2], -v.v[3] },
            };
        }

        /// Element-wise multiplication.
        pub inline fn mul(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = v.v * u.v,
            };
        }

        /// Scalar multiplication.
        pub inline fn mulScalar(v: *const Vector4(T), s: T) Vector4(T) {
            return Vector4(T){
                .v = v.v * @Vector(4, T){ s, s, s, s },
            };
        }

        /// Element-wise division.
        pub inline fn div(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = v.v / u.v,
            };
        }

        /// Scalar division.
        pub inline fn divScalar(v: *const Vector4(T), s: T) Vector4(T) {
            return Vector4(T){
                .v = v.v / @Vector(4, T){ s, s, s, s },
            };
        }

        /// Computes the modulus of a vector.
        pub inline fn mod(v: *const Vector4(T)) T {
            return @sqrt(v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2] + v.v[3] * v.v[3]);
        }

        /// Computes the squared modulus of a vector.
        pub inline fn mod2(v: *const Vector4(T)) T {
            return v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2] + v.v[3] * v.v[3];
        }

        /// Normalizes a vector.
        pub inline fn normalize(v: *const Vector4(T)) Vector4(T) {
            return v.div(&splat(v.mod()));
        }

        /// Computes the dot product between two vectors.
        pub inline fn dot(v: *const Vector4(T), u: *const Vector4(T)) T {
            return v.v[0] * u.v[0] + v.v[1] * u.v[1] + v.v[2] * u.v[2] + v.v[3] * u.v[3];
        }

        /// Computes the distance between two vectors.
        pub inline fn distance(v: *const Vector4(T), u: *const Vector4(T)) T {
            return v.sub(u).mod();
        }

        /// Computes the squared distance between two vectors.
        pub inline fn distance2(v: *const Vector4(T), u: *const Vector4(T)) T {
            return v.sub(u).mod2();
        }

        /// Linearly interpolates two vectors.
        pub inline fn lerp(v: *const Vector4(T), u: *const Vector4(T), t: T) Vector4(T) {
            return v.mulScalar(1 - t).add(&u.mulScalar(t));
        }

        /// Clamps element-wise a vector between a maximum and a minimum.
        pub inline fn clamp(v: *const Vector4(T), mn: T, mx: T) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ @min(@max(v.v[0], mn), mx), @min(@max(v.v[1], mn), mx), @min(@max(v.v[2], mn), mx), @min(@max(v.v[3], mn), mx) },
            };
        }

        /// Computes the angle between two vectors.
        pub inline fn angle(v: *const Vector4(T), u: *const Vector4(T)) T {
            return std.math.acos((v.dot(u)) / (v.mod() * u.mod()));
        }

        /// Computes the projection of the first vector onto the second.
        pub inline fn project(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return u.mul(&splat(v.dot(u) / u.mod2()));
        }

        /// Reflects a vector against a normal.
        pub inline fn reflect(v: *const Vector4(T), normal: *const Vector4(T)) Vector4(T) {
            return v.sub(&normal.mul(&splat(2 * v.dot(normal) / normal.mod2())));
        }

        /// Refracts a vector against a normal.
        pub inline fn refract(v: *const Vector4(T), normal: *const Vector4(T), eta: T) Vector4(T) {
            const v_norm = v.normalize();
            const k = 1 - eta * eta * (1 - v_norm.dot(normal) * v_norm.dot(normal));

            if (k < 0) {
                return splat(0);
            }

            return v_norm.mulScalar(eta).sub(&normal.mulScalar(eta * v_norm.dot(normal) + @sqrt(k)));
        }

        /// Returns a vector with the maximum values between two vectors.
        pub inline fn max(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ @max(v.v[0], u.v[0]), @max(v.v[1], u.v[1]), @max(v.v[2], u.v[2]), @max(v.v[3], u.v[3]) },
            };
        }

        /// Returns a vector with the minimum values between two vectors.
        pub inline fn min(v: *const Vector4(T), u: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ @min(v.v[0], u.v[0]), @min(v.v[1], u.v[1]), @min(v.v[2], u.v[2]), @min(v.v[3], u.v[3]) },
            };
        }

        /// Returns a vector with the absolute values of the input vector.
        pub inline fn abs(v: *const Vector4(T)) Vector4(T) {
            return Vector4(T){
                .v = @Vector(4, T){ @abs(v.v[0]), @abs(v.v[1]), @abs(v.v[2]), @abs(v.v[3]) },
            };
        }

        /// Compares two vectors and returns true if they are equal.
        pub inline fn equal(v: *const Vector4(T), u: *const Vector4(T)) bool {
            return v.v[0] == u.v[0] and v.v[1] == u.v[1] and v.v[2] == u.v[2] and v.v[3] == u.v[3];
        }

        /// Compares two vectors and returns true if they are approximately equal.
        pub inline fn approxEqual(v: *const Vector4(T), u: *const Vector4(T), tolerance: T) bool {
            return @abs(v.v[0] - u.v[0]) <= tolerance and @abs(v.v[1] - u.v[1]) <= tolerance and @abs(v.v[2] - u.v[2]) <= tolerance and @abs(v.v[3] - u.v[3]) <= tolerance;
        }
    };
}

test "add_vector2" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector2(f32).init(3, 4);
    const r = v.add(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(4, 6), 0.0001));
}

test "addScalar_vector2" {
    const v = Vector2(f32).init(1, 2);
    const r = v.addScalar(3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(4, 5), 0.0001));
}

test "sub_vector2" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector2(f32).init(3, 4);
    const r = v.sub(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(-2, -2), 0.0001));
}

test "subScalar_vector2" {
    const v = Vector2(f32).init(1, 2);
    const r = v.subScalar(3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(-2, -1), 0.0001));
}

test "negate_vector2" {
    const v = Vector2(f32).init(1, 2);
    const r = v.negate();

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(-1, -2), 0.0001));
}

test "mul_vector2" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector2(f32).init(3, 4);
    const r = v.mul(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(3, 8), 0.0001));
}

test "mulScalar_vector2" {
    const v = Vector2(f32).init(1, 2);
    const r = v.mulScalar(3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(3, 6), 0.0001));
}

test "div_vector2" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector2(f32).init(3, 4);
    const r = v.div(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(0.3333, 0.5), 0.0001));
}

test "divScalar_vector2" {
    const v = Vector2(f32).init(1, 2);
    const r = v.divScalar(3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(0.3333, 0.6666), 0.0001));
}

test "mod_vector2" {
    const v = Vector2(f32).init(3, 4);
    const r = v.mod();

    try std.testing.expectApproxEqRel(5, r, 0.001);
}

test "mod2_vector2" {
    const v = Vector2(f32).init(3, 4);
    const r = v.mod2();

    try std.testing.expectApproxEqRel(25, r, 0.001);
}

test "normalize_vector2" {
    const v = Vector2(f32).init(3, 4);
    const r = v.normalize();

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(0.6, 0.8), 0.0001));
}

test "dot_vector2" {
    const v = Vector2(f32).init(3, 4);
    const w = Vector2(f32).init(5, 6);
    const r = v.dot(&w);

    try std.testing.expectApproxEqRel(39, r, 0.001);
}

test "distance_vector2" {
    const v = Vector2(f32).init(3, 4);
    const w = Vector2(f32).init(5, 6);
    const r = v.distance(&w);

    try std.testing.expectApproxEqRel(2.8284, r, 0.001);
}

test "distance2_vector2" {
    const v = Vector2(f32).init(3, 4);
    const w = Vector2(f32).init(5, 6);
    const r = v.distance2(&w);

    try std.testing.expectApproxEqRel(8, r, 0.001);
}

test "lerp_vector2" {
    const v = Vector2(f32).init(3, 4);
    const w = Vector2(f32).init(5, 6);
    const r = v.lerp(&w, 0.5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(4, 5), 0.0001));
}

test "clamp_vector2" {
    const v = Vector2(f32).init(3, 4);
    const r = v.clamp(2, 3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(3, 3), 0.0001));
}

test "angle_vector2" {
    const v = Vector2(f32).init(3, 4);
    const w = Vector2(f32).init(5, 6);
    const r = v.angle(&w);

    try std.testing.expectApproxEqRel(0.0512, r, 0.001);
}

test "project_vector2" {
    const v = Vector2(f32).init(3, 4);
    const w = Vector2(f32).init(5, 6);
    const r = v.project(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(3.1967, 3.8361), 0.0001));
}

test "reflect_vector2" {
    const v = Vector2(f32).init(3, 4);
    const normal = Vector2(f32).init(1, 0);
    const r = v.reflect(&normal);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(-3, 4), 0.0001));
}

test "refract_vector2" {
    const v = Vector2(f32).init(3, 4);
    const normal = Vector2(f32).init(1, 0);
    const r = v.refract(&normal, 1.1);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(-0.4749, 0.88), 0.0001));
}

test "max_vector2" {
    const v = Vector2(f32).init(3, 6);
    const w = Vector2(f32).init(5, 4);
    const r = v.max(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(5, 6), 0.0001));
}

test "min_vector2" {
    const v = Vector2(f32).init(3, 6);
    const w = Vector2(f32).init(5, 4);
    const r = v.min(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(3, 4), 0.0001));
}

test "abs_vector2" {
    const v = Vector2(f32).init(-3, -6);
    const r = v.abs();

    try std.testing.expectEqual(true, r.approxEqual(&Vector2(f32).init(3, 6), 0.0001));
}

test "add_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.add(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(5, 7, 9), 0.0001));
}

test "addScalar_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.addScalar(4);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(5, 6, 7), 0.0001));
}

test "sub_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.sub(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(-3, -3, -3), 0.0001));
}

test "subScalar_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.subScalar(4);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(-3, -2, -1), 0.0001));
}

test "negate_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.negate();

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(-1, -2, -3), 0.0001));
}

test "mul_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.mul(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(4, 10, 18), 0.0001));
}

test "mulScalar_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.mulScalar(4);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(4, 8, 12), 0.0001));
}

test "div_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.div(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(0.25, 0.4, 0.5), 0.0001));
}

test "divScalar_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.divScalar(4);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(0.25, 0.5, 0.75), 0.0001));
}

test "mod_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.mod();

    try std.testing.expectApproxEqRel(3.7416, r, 0.001);
}

test "mod2_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.mod2();

    try std.testing.expectApproxEqRel(14, r, 0.001);
}

test "normalize_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.normalize();

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(0.2672, 0.5345, 0.8017), 0.0001));
}

test "dot_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.dot(&w);

    try std.testing.expectApproxEqRel(32, r, 0.001);
}

test "cross_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.cross(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(-3, 6, -3), 0.0001));
}

test "distance_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.distance(&w);

    try std.testing.expectApproxEqRel(5.1961, r, 0.001);
}

test "distance2_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.distance2(&w);

    try std.testing.expectApproxEqRel(27, r, 0.001);
}

test "lerp_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.lerp(&w, 0.5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(2.5, 3.5, 4.5), 0.0001));
}

test "clamp_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const r = v.clamp(2, 3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(2, 2, 3), 0.0001));
}

test "angle_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.angle(&w);

    try std.testing.expectApproxEqRel(0.2257, r, 0.001);
}

test "project_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector3(f32).init(4, 5, 6);
    const r = v.project(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(1.6623, 2.0779, 2.4935), 0.0001));
}

test "reflect_vector3" {
    const v = Vector3(f32).init(1, 2, 3);
    const normal = Vector3(f32).init(1, 0, 0);
    const r = v.reflect(&normal);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(-1, 2, 3), 0.0001));
}

test "refract_vector3" {
    const v = Vector3(f32).init(3, 2, 1);
    const normal = Vector3(f32).init(1, 0, 1);
    const r = v.refract(&normal, 1.1);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(-1.3769, 0.5879, -1.9649), 0.0001));
}

test "max_vector3" {
    const v = Vector3(f32).init(1, 5, 3);
    const w = Vector3(f32).init(4, 2, 6);
    const r = v.max(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(4, 5, 6), 0.0001));
}

test "min_vector3" {
    const v = Vector3(f32).init(1, 5, 3);
    const w = Vector3(f32).init(4, 2, 6);
    const r = v.min(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(1, 2, 3), 0.0001));
}

test "abs_vector3" {
    const v = Vector3(f32).init(-1, -5, -3);
    const r = v.abs();

    try std.testing.expectEqual(true, r.approxEqual(&Vector3(f32).init(1, 5, 3), 0.0001));
}

test "add_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.add(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(6, 8, 10, 12), 0.0001));
}

test "addScalar_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.addScalar(5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(6, 7, 8, 9), 0.0001));
}

test "sub_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.sub(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(-4, -4, -4, -4), 0.0001));
}

test "subScalar_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.subScalar(5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(-4, -3, -2, -1), 0.0001));
}

test "negate_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.negate();

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(-1, -2, -3, -4), 0.0001));
}

test "mul_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.mul(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(5, 12, 21, 32), 0.0001));
}

test "mulScalar_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.mulScalar(5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(5, 10, 15, 20), 0.0001));
}

test "div_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.div(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(0.2, 0.3333, 0.4285, 0.5), 0.0001));
}

test "divScalar_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.divScalar(5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(0.2, 0.4, 0.6, 0.8), 0.0001));
}

test "mod_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.mod();

    try std.testing.expectApproxEqRel(5.4772, r, 0.001);
}

test "mod2_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.mod2();

    try std.testing.expectApproxEqRel(30, r, 0.001);
}

test "normalize_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.normalize();

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(0.1826, 0.3651, 0.5477, 0.7302), 0.0001));
}

test "dot_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.dot(&w);

    try std.testing.expectApproxEqRel(70, r, 0.001);
}

test "distance_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.distance(&w);

    try std.testing.expectApproxEqRel(8, r, 0.001);
}

test "distance2_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.distance2(&w);

    try std.testing.expectApproxEqRel(64, r, 0.001);
}

test "lerp_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.lerp(&w, 0.5);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(3, 4, 5, 6), 0.0001));
}

test "clamp_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const r = v.clamp(2, 3);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(2, 2, 3, 3), 0.0001));
}

test "angle_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.angle(&w);

    try std.testing.expectApproxEqRel(0.2501, r, 0.001);
}

test "project_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const w = Vector4(f32).init(5, 6, 7, 8);
    const r = v.project(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(2.0114, 2.4137, 2.8160, 3.2183), 0.0001));
}

test "reflect_vector4" {
    const v = Vector4(f32).init(1, 2, 3, 4);
    const normal = Vector4(f32).init(1, 0, 0, 0);
    const r = v.reflect(&normal);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(-1, 2, 3, 4), 0.0001));
}

test "refract_vector4" {
    const v = Vector4(f32).init(3, 2, 1, 4);
    const normal = Vector4(f32).init(1, 0, 1, 0);
    const r = v.refract(&normal, 1.1);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(-0.8606, 0.4016, -1.2622, 0.80332), 0.0001));
}

test "max_vector4" {
    const v = Vector4(f32).init(1, 5, 3, 4);
    const w = Vector4(f32).init(4, 2, 6, 8);
    const r = v.max(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(4, 5, 6, 8), 0.0001));
}

test "min_vector4" {
    const v = Vector4(f32).init(1, 5, 3, 4);
    const w = Vector4(f32).init(4, 2, 6, 8);
    const r = v.min(&w);

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(1, 2, 3, 4), 0.0001));
}

test "abs_vector4" {
    const v = Vector4(f32).init(-1, -5, -3, -4);
    const r = v.abs();

    try std.testing.expectEqual(true, r.approxEqual(&Vector4(f32).init(1, 5, 3, 4), 0.0001));
}
