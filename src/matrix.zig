const std = @import("std");
const zgm = @import("zgm.zig");
const Vector2 = zgm.Vector2;
const Vector3 = zgm.Vector3;
const Vector4 = zgm.Vector4;

/// A 2x2 matrix.
pub fn Matrix2x2(comptime T: type) type {
    return struct {
        m: [2]Vector2(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m10: T, m11: T) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(m00, m10),
                Vector2(T).init(m01, m11),
            } };
        }

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector2(T), r1: *const Vector2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(r0.v[0], r1.v[0]),
                Vector2(T).init(r0.v[1], r1.v[1]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initColumn(c0: *const Vector2(T), c1: *const Vector2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(c0.v[0], c0.v[1]),
                Vector2(T).init(c1.v[0], c1.v[1]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).splat(s),
                Vector2(T).splat(s),
            } };
        }

        /// Identity matrix.
        pub const eye = Matrix2x2(T).init(1, 0, 0, 1);

        /// Zero matrix.
        pub const zero = Matrix2x2(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix2x2(T), r: usize, c: usize) T {
            return A.m[c].v[r];
        }

        /// Get the row vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        ///
        /// **Returns**:
        /// - The row vector at the given index.
        pub inline fn row(A: *const Matrix2x2(T), r: usize) Vector2(T) {
            return Vector2(T).init(A.m[0].v[r], A.m[1].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix2x2(T), c: usize) Vector2(T) {
            return Vector2(T).init(A.m[c].v[0], A.m[c].v[1]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrix.
        pub inline fn add(A: *const Matrix2x2(T), B: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).add(&A.m[0], &B.m[0]),
                Vector2(T).add(&A.m[1], &B.m[1]),
            } };
        }

        /// Adds a scalar to all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new matrix with the scalar added to all components.
        pub inline fn addScalar(A: *const Matrix2x2(T), s: T) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).addScalar(&A.m[0], s),
                Vector2(T).addScalar(&A.m[1], s),
            } };
        }

        /// Element-wise subtraction.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise difference of the input
        ///   matrices.
        pub inline fn sub(A: *const Matrix2x2(T), B: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).sub(&A.m[0], &B.m[0]),
                Vector2(T).sub(&A.m[1], &B.m[1]),
            } };
        }

        /// Subtracts a scalar from all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new matrix with the scalar subtracted from all components.
        pub inline fn subScalar(A: *const Matrix2x2(T), s: T) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).subScalar(&A.m[0], s),
                Vector2(T).subScalar(&A.m[1], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).negate(&A.m[0]),
                Vector2(T).negate(&A.m[1]),
            } };
        }

        /// Element-wise multiplication.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise product of the input matrices.
        pub inline fn mul(A: *const Matrix2x2(T), B: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).mul(&A.m[0], &B.m[0]),
                Vector2(T).mul(&A.m[1], &B.m[1]),
            } };
        }

        /// Multiplies all components of a matrix by a scalar.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new matrix with the components multiplied by the scalar.
        pub inline fn mulScalar(A: *const Matrix2x2(T), s: T) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).mulScalar(&A.m[0], s),
                Vector2(T).mulScalar(&A.m[1], s),
            } };
        }

        /// Matrix2x2-Vector2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector2(A: *const Matrix2x2(T), v: Vector2(T)) Vector2(T) {
            return Vector2(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1],
            );
        }

        /// Matrix2x2-Matrix2x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x2(A: *const Matrix2x2(T), B: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                ),
            } };
        }

        /// Matrix2x2-Matrix2x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x3(A: *const Matrix2x2(T), B: *const Matrix2x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1],
                ),
            } };
        }

        /// Matrix2x2-Matrix2x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x4(A: *const Matrix2x2(T), B: *const Matrix2x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1],
                ),
            } };
        }

        /// Element-wise division.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise division of the input matrices.
        pub inline fn div(A: *const Matrix2x2(T), B: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).div(&A.m[0], &B.m[0]),
                Vector2(T).div(&A.m[1], &B.m[1]),
            } };
        }

        /// Divides all components of a matrix by a scalar.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `s`: The scalar.
        ///
        /// **Returns**:
        /// - A new matrix with the components divided by the scalar.
        pub inline fn divScalar(A: *const Matrix2x2(T), s: T) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).divScalar(&A.m[0], s),
                Vector2(T).divScalar(&A.m[1], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(A.m[0].v[0], A.m[0].v[1]),
                Vector2(T).init(A.m[1].v[0], A.m[1].v[1]),
            } };
        }

        /// Determinant of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The determinant of the matrix.
        pub inline fn det(A: *const Matrix2x2(T)) T {
            return A.m[0].v[0] * A.m[1].v[1] - A.m[0].v[1] * A.m[1].v[0];
        }

        /// Inverse of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The inverse of the matrix.
        /// - If the matrix is not invertible, an error is returned.
        pub inline fn inverse(A: *const Matrix2x2(T)) !Matrix2x2(T) {
            const determinant = Matrix2x2(T).det(A);
            if (determinant == 0) {
                return error.SingularMatrix;
            }

            const invDet = 1 / determinant;
            return Matrix2x2(T).init(
                A.m[1].v[1] * invDet,
                -A.m[1].v[0] * invDet,
                -A.m[0].v[1] * invDet,
                A.m[0].v[0] * invDet,
            );
        }
    };
}

/// A 2x3 matrix.
pub fn Matrix2x3(comptime T: type) type {
    return struct {
        m: [3]Vector2(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m02`: The value at row 0, column 2.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m12`: The value at row 1, column 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m02: T, m10: T, m11: T, m12: T) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(m00, m10),
                Vector2(T).init(m01, m11),
                Vector2(T).init(m02, m12),
            } };
        }
    };
}

/// A 2x4 matrix.
pub fn Matrix2x4(comptime T: type) type {
    return struct {
        m: [4]Vector2(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m02`: The value at row 0, column 2.
        /// - `m03`: The value at row 0, column 3.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m12`: The value at row 1, column 2.
        /// - `m13`: The value at row 1, column 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m02: T, m03: T, m10: T, m11: T, m12: T, m13: T) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(m00, m10),
                Vector2(T).init(m01, m11),
                Vector2(T).init(m02, m12),
                Vector2(T).init(m03, m13),
            } };
        }
    };
}

/// A 3x2 matrix.
pub fn Matrix3x2(comptime T: type) type {
    return struct {
        m: [2]Vector3(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m20`: The value at row 2, column 0.
        /// - `m21`: The value at row 2, column 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m10: T, m11: T, m20: T, m21: T) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(m00, m10, m20),
                Vector3(T).init(m01, m11, m21),
            } };
        }
    };
}

/// A 3x3 matrix.
pub fn Matrix3x3(comptime T: type) type {
    return struct {
        m: [3]Vector3(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m02`: The value at row 0, column 2.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m12`: The value at row 1, column 2.
        /// - `m20`: The value at row 2, column 0.
        /// - `m21`: The value at row 2, column 1.
        /// - `m22`: The value at row 2, column 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m02: T, m10: T, m11: T, m12: T, m20: T, m21: T, m22: T) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(m00, m10, m20),
                Vector3(T).init(m01, m11, m21),
                Vector3(T).init(m02, m12, m22),
            } };
        }
    };
}

/// A 3x4 matrix.
pub fn Matrix3x4(comptime T: type) type {
    return struct {
        m: [4]Vector3(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m02`: The value at row 0, column 2.
        /// - `m03`: The value at row 0, column 3.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m12`: The value at row 1, column 2.
        /// - `m13`: The value at row 1, column 3.
        /// - `m20`: The value at row 2, column 0.
        /// - `m21`: The value at row 2, column 1.
        /// - `m22`: The value at row 2, column 2.
        /// - `m23`: The value at row 2, column 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m02: T, m03: T, m10: T, m11: T, m12: T, m13: T, m20: T, m21: T, m22: T, m23: T) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).init(m00, m10, m20),
                Vector3(T).init(m01, m11, m21),
                Vector3(T).init(m02, m12, m22),
                Vector3(T).init(m03, m13, m23),
            } };
        }
    };
}

/// A 4x2 matrix.
pub fn Matrix4x2(comptime T: type) type {
    return struct {
        m: [2]Vector4(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m20`: The value at row 2, column 0.
        /// - `m21`: The value at row 2, column 1.
        /// - `m30`: The value at row 3, column 0.
        /// - `m31`: The value at row 3, column 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m10: T, m11: T, m20: T, m21: T, m30: T, m31: T) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(m00, m10, m20, m30),
                Vector4(T).init(m01, m11, m21, m31),
            } };
        }
    };
}

/// A 4x3 matrix.
pub fn Matrix4x3(comptime T: type) type {
    return struct {
        m: [3]Vector4(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m02`: The value at row 0, column 2.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m12`: The value at row 1, column 2.
        /// - `m20`: The value at row 2, column 0.
        /// - `m21`: The value at row 2, column 1.
        /// - `m22`: The value at row 2, column 2.
        /// - `m30`: The value at row 3, column 0.
        /// - `m31`: The value at row 3, column 1.
        /// - `m32`: The value at row 3, column 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m02: T, m10: T, m11: T, m12: T, m20: T, m21: T, m22: T, m30: T, m31: T, m32: T) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(m00, m10, m20, m30),
                Vector4(T).init(m01, m11, m21, m31),
                Vector4(T).init(m02, m12, m22, m32),
            } };
        }
    };
}

/// A 4x4 matrix.
pub fn Matrix4x4(comptime T: type) type {
    return struct {
        m: [4]Vector4(T),

        /// Initialize a matrix with the given values.
        ///
        /// **Parameters**:
        /// - `m00`: The value at row 0, column 0.
        /// - `m01`: The value at row 0, column 1.
        /// - `m02`: The value at row 0, column 2.
        /// - `m03`: The value at row 0, column 3.
        /// - `m10`: The value at row 1, column 0.
        /// - `m11`: The value at row 1, column 1.
        /// - `m12`: The value at row 1, column 2.
        /// - `m13`: The value at row 1, column 3.
        /// - `m20`: The value at row 2, column 0.
        /// - `m21`: The value at row 2, column 1.
        /// - `m22`: The value at row 2, column 2.
        /// - `m23`: The value at row 2, column 3.
        /// - `m30`: The value at row 3, column 0.
        /// - `m31`: The value at row 3, column 1.
        /// - `m32`: The value at row 3, column 2.
        /// - `m33`: The value at row 3, column 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given values.
        pub inline fn init(m00: T, m01: T, m02: T, m03: T, m10: T, m11: T, m12: T, m13: T, m20: T, m21: T, m22: T, m23: T, m30: T, m31: T, m32: T, m33: T) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(m00, m10, m20, m30),
                Vector4(T).init(m01, m11, m21, m31),
                Vector4(T).init(m02, m12, m22, m32),
                Vector4(T).init(m03, m13, m23, m33),
            } };
        }
    };
}
