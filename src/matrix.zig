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
        pub inline fn initCol(c0: *const Vector2(T), c1: *const Vector2(T)) Matrix2x2(T) {
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
        pub const identity = Matrix2x2(T).init(1, 0, 0, 1);

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
        pub inline fn mulVector2(A: *const Matrix2x2(T), v: *const Vector2(T)) Vector2(T) {
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
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix2x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(A.m[0].v[0], A.m[1].v[0]),
                Vector2(T).init(A.m[0].v[1], A.m[1].v[1]),
            } };
        }

        /// Determinant of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The determinant of the matrix.
        pub inline fn determinant(A: *const Matrix2x2(T)) T {
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
            const det = Matrix2x2(T).determinant(A);
            if (det == 0) {
                return error.SingularMatrix;
            }

            const invDet = 1 / det;
            return Matrix2x2(T).init(
                A.m[1].v[1] * invDet,
                -A.m[1].v[0] * invDet,
                -A.m[0].v[1] * invDet,
                A.m[0].v[0] * invDet,
            );
        }

        /// Trace of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The trace of the matrix.
        pub inline fn trace(A: *const Matrix2x2(T)) T {
            return A.m[0].v[0] + A.m[1].v[1];
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix2x2(T), B: *const Matrix2x2(T)) bool {
            return Vector2(T).equal(&A.m[0], &B.m[0]) and Vector2(T).equal(&A.m[1], &B.m[1]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix2x2(T), B: *const Matrix2x2(T), tolerance: T) bool {
            return Vector2(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector2(T).approxEqual(&A.m[1], &B.m[1], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector3(T), r1: *const Vector3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(r0.v[0], r1.v[0]),
                Vector2(T).init(r0.v[1], r1.v[1]),
                Vector2(T).init(r0.v[2], r1.v[2]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        /// - `c2`: The vector for column 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initCol(c0: *const Vector2(T), c1: *const Vector2(T), c2: *const Vector2(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(c0.v[0], c0.v[1]),
                Vector2(T).init(c1.v[0], c1.v[1]),
                Vector2(T).init(c2.v[0], c2.v[1]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).splat(s),
                Vector2(T).splat(s),
                Vector2(T).splat(s),
            } };
        }

        /// Zero matrix.
        pub const zero = Matrix2x3(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix2x3(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix2x3(T), r: usize) Vector3(T) {
            return Vector3(T).init(A.m[0].v[r], A.m[1].v[r], A.m[2].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix2x3(T), c: usize) Vector2(T) {
            return Vector2(T).init(A.m[c].v[0], A.m[c].v[1]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix2x3(T), B: *const Matrix2x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).add(&A.m[0], &B.m[0]),
                Vector2(T).add(&A.m[1], &B.m[1]),
                Vector2(T).add(&A.m[2], &B.m[2]),
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
        pub inline fn addScalar(A: *const Matrix2x3(T), s: T) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).addScalar(&A.m[0], s),
                Vector2(T).addScalar(&A.m[1], s),
                Vector2(T).addScalar(&A.m[2], s),
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
        pub inline fn sub(A: *const Matrix2x3(T), B: *const Matrix2x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).sub(&A.m[0], &B.m[0]),
                Vector2(T).sub(&A.m[1], &B.m[1]),
                Vector2(T).sub(&A.m[2], &B.m[2]),
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
        pub inline fn subScalar(A: *const Matrix2x3(T), s: T) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).subScalar(&A.m[0], s),
                Vector2(T).subScalar(&A.m[1], s),
                Vector2(T).subScalar(&A.m[2], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix2x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).negate(&A.m[0]),
                Vector2(T).negate(&A.m[1]),
                Vector2(T).negate(&A.m[2]),
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
        pub inline fn mul(A: *const Matrix2x3(T), B: *const Matrix2x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).mul(&A.m[0], &B.m[0]),
                Vector2(T).mul(&A.m[1], &B.m[1]),
                Vector2(T).mul(&A.m[2], &B.m[2]),
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
        pub inline fn mulScalar(A: *const Matrix2x3(T), s: T) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).mulScalar(&A.m[0], s),
                Vector2(T).mulScalar(&A.m[1], s),
                Vector2(T).mulScalar(&A.m[2], s),
            } };
        }

        /// Matrix2x3-Vector3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector3(A: *const Matrix2x3(T), v: *const Vector3(T)) Vector2(T) {
            return Vector2(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2],
            );
        }

        /// Matrix2x3-Matrix3x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x2(A: *const Matrix2x3(T), B: *const Matrix3x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                ),
            } };
        }

        /// Matrix2x3-Matrix3x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x3(A: *const Matrix2x3(T), B: *const Matrix3x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2],
                ),
            } };
        }

        /// Matrix2x3-Matrix3x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x4(A: *const Matrix2x3(T), B: *const Matrix3x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[2].v[0] * B.m[3].v[2],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[2].v[1] * B.m[3].v[2],
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
        pub inline fn div(A: *const Matrix2x3(T), B: *const Matrix2x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).div(&A.m[0], &B.m[0]),
                Vector2(T).div(&A.m[1], &B.m[1]),
                Vector2(T).div(&A.m[2], &B.m[2]),
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
        pub inline fn divScalar(A: *const Matrix2x3(T), s: T) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).divScalar(&A.m[0], s),
                Vector2(T).divScalar(&A.m[1], s),
                Vector2(T).divScalar(&A.m[2], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix2x3(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(A.m[0].v[0], A.m[1].v[0], A.m[2].v[0]),
                Vector3(T).init(A.m[0].v[1], A.m[1].v[1], A.m[2].v[1]),
            } };
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix2x3(T), B: *const Matrix2x3(T)) bool {
            return Vector2(T).equal(&A.m[0], &B.m[0]) and Vector2(T).equal(&A.m[1], &B.m[1]) and Vector2(T).equal(&A.m[2], &B.m[2]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix2x3(T), B: *const Matrix2x3(T), tolerance: T) bool {
            return Vector2(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector2(T).approxEqual(&A.m[1], &B.m[1], tolerance) and Vector2(T).approxEqual(&A.m[2], &B.m[2], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector4(T), r1: *const Vector4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(r0.v[0], r1.v[0]),
                Vector2(T).init(r0.v[1], r1.v[1]),
                Vector2(T).init(r0.v[2], r1.v[2]),
                Vector2(T).init(r0.v[3], r1.v[3]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        /// - `c2`: The vector for column 2.
        /// - `c3`: The vector for column 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initCol(c0: *const Vector2(T), c1: *const Vector2(T), c2: *const Vector2(T), c3: *const Vector2(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(c0.v[0], c0.v[1]),
                Vector2(T).init(c1.v[0], c1.v[1]),
                Vector2(T).init(c2.v[0], c2.v[1]),
                Vector2(T).init(c3.v[0], c3.v[1]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).splat(s),
                Vector2(T).splat(s),
                Vector2(T).splat(s),
                Vector2(T).splat(s),
            } };
        }

        /// Zero matrix.
        pub const zero = Matrix2x4(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix2x4(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix2x4(T), r: usize) Vector4(T) {
            return Vector4(T).init(A.m[0].v[r], A.m[1].v[r], A.m[2].v[r], A.m[3].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix2x4(T), c: usize) Vector2(T) {
            return Vector2(T).init(A.m[c].v[0], A.m[c].v[1]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix2x4(T), B: *const Matrix2x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).add(&A.m[0], &B.m[0]),
                Vector2(T).add(&A.m[1], &B.m[1]),
                Vector2(T).add(&A.m[2], &B.m[2]),
                Vector2(T).add(&A.m[3], &B.m[3]),
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
        pub inline fn addScalar(A: *const Matrix2x4(T), s: T) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).addScalar(&A.m[0], s),
                Vector2(T).addScalar(&A.m[1], s),
                Vector2(T).addScalar(&A.m[2], s),
                Vector2(T).addScalar(&A.m[3], s),
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
        pub inline fn sub(A: *const Matrix2x4(T), B: *const Matrix2x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).sub(&A.m[0], &B.m[0]),
                Vector2(T).sub(&A.m[1], &B.m[1]),
                Vector2(T).sub(&A.m[2], &B.m[2]),
                Vector2(T).sub(&A.m[3], &B.m[3]),
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
        pub inline fn subScalar(A: *const Matrix2x4(T), s: T) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).subScalar(&A.m[0], s),
                Vector2(T).subScalar(&A.m[1], s),
                Vector2(T).subScalar(&A.m[2], s),
                Vector2(T).subScalar(&A.m[3], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix2x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).negate(&A.m[0]),
                Vector2(T).negate(&A.m[1]),
                Vector2(T).negate(&A.m[2]),
                Vector2(T).negate(&A.m[3]),
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
        pub inline fn mul(A: *const Matrix2x4(T), B: *const Matrix2x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).mul(&A.m[0], &B.m[0]),
                Vector2(T).mul(&A.m[1], &B.m[1]),
                Vector2(T).mul(&A.m[2], &B.m[2]),
                Vector2(T).mul(&A.m[3], &B.m[3]),
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
        pub inline fn mulScalar(A: *const Matrix2x4(T), s: T) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).mulScalar(&A.m[0], s),
                Vector2(T).mulScalar(&A.m[1], s),
                Vector2(T).mulScalar(&A.m[2], s),
                Vector2(T).mulScalar(&A.m[3], s),
            } };
        }

        /// Matrix2x4-Vector4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector4(A: *const Matrix2x4(T), v: *const Vector4(T)) Vector2(T) {
            return Vector2(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2] + A.m[3].v[0] * v.v[3],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2] + A.m[3].v[1] * v.v[3],
            );
        }

        /// Matrix2x4-Matrix4x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x2(A: *const Matrix2x4(T), B: *const Matrix4x2(T)) Matrix2x2(T) {
            return Matrix2x2(T){ .m = [2]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                ),
            } };
        }

        /// Matrix2x4-Matrix4x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x3(A: *const Matrix2x4(T), B: *const Matrix4x3(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2] + A.m[3].v[0] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2] + A.m[3].v[1] * B.m[2].v[3],
                ),
            } };
        }

        /// Matrix2x4-Matrix4x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x4(A: *const Matrix2x4(T), B: *const Matrix4x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2] + A.m[3].v[0] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2] + A.m[3].v[1] * B.m[2].v[3],
                ),
                Vector2(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[2].v[0] * B.m[3].v[2] + A.m[3].v[0] * B.m[3].v[3],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[2].v[1] * B.m[3].v[2] + A.m[3].v[1] * B.m[3].v[3],
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
        pub inline fn div(A: *const Matrix2x4(T), B: *const Matrix2x4(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).div(&A.m[0], &B.m[0]),
                Vector2(T).div(&A.m[1], &B.m[1]),
                Vector2(T).div(&A.m[2], &B.m[2]),
                Vector2(T).div(&A.m[3], &B.m[3]),
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
        pub inline fn divScalar(A: *const Matrix2x4(T), s: T) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).divScalar(&A.m[0], s),
                Vector2(T).divScalar(&A.m[1], s),
                Vector2(T).divScalar(&A.m[2], s),
                Vector2(T).divScalar(&A.m[3], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix2x4(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(A.m[0].v[0], A.m[1].v[0], A.m[2].v[0], A.m[3].v[0]),
                Vector4(T).init(A.m[0].v[1], A.m[1].v[1], A.m[2].v[1], A.m[3].v[1]),
            } };
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix2x4(T), B: *const Matrix2x4(T)) bool {
            return Vector2(T).equal(&A.m[0], &B.m[0]) and Vector2(T).equal(&A.m[1], &B.m[1]) and Vector2(T).equal(&A.m[2], &B.m[2]) and Vector2(T).equal(&A.m[3], &B.m[3]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix2x4(T), B: *const Matrix2x4(T), tolerance: T) bool {
            return Vector2(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector2(T).approxEqual(&A.m[1], &B.m[1], tolerance) and Vector2(T).approxEqual(&A.m[2], &B.m[2], tolerance) and Vector2(T).approxEqual(&A.m[3], &B.m[3], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        /// - `r2`: The vector for row 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector2(T), r1: *const Vector2(T), r2: *const Vector2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(r0.v[0], r1.v[0], r2.v[0]),
                Vector3(T).init(r0.v[1], r1.v[1], r2.v[1]),
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
        pub inline fn initCol(c0: *const Vector3(T), c1: *const Vector3(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(c0.v[0], c0.v[1], c0.v[2]),
                Vector3(T).init(c1.v[0], c1.v[1], c1.v[2]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).splat(s),
                Vector3(T).splat(s),
            } };
        }

        /// Zero matrix.
        pub const zero = Matrix3x2(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix3x2(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix3x2(T), r: usize) Vector2(T) {
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
        pub inline fn col(A: *const Matrix3x2(T), c: usize) Vector3(T) {
            return Vector3(T).init(A.m[c].v[0], A.m[c].v[1], A.m[c].v[2]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix3x2(T), B: *const Matrix3x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).add(&A.m[0], &B.m[0]),
                Vector3(T).add(&A.m[1], &B.m[1]),
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
        pub inline fn addScalar(A: *const Matrix3x2(T), s: T) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).addScalar(&A.m[0], s),
                Vector3(T).addScalar(&A.m[1], s),
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
        pub inline fn sub(A: *const Matrix3x2(T), B: *const Matrix3x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).sub(&A.m[0], &B.m[0]),
                Vector3(T).sub(&A.m[1], &B.m[1]),
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
        pub inline fn subScalar(A: *const Matrix3x2(T), s: T) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).subScalar(&A.m[0], s),
                Vector3(T).subScalar(&A.m[1], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix3x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).negate(&A.m[0]),
                Vector3(T).negate(&A.m[1]),
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
        pub inline fn mul(A: *const Matrix3x2(T), B: *const Matrix3x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).mul(&A.m[0], &B.m[0]),
                Vector3(T).mul(&A.m[1], &B.m[1]),
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
        pub inline fn mulScalar(A: *const Matrix3x2(T), s: T) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).mulScalar(&A.m[0], s),
                Vector3(T).mulScalar(&A.m[1], s),
            } };
        }

        /// Matrix3x2-Vector2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector2(A: *const Matrix3x2(T), v: *const Vector2(T)) Vector3(T) {
            return Vector3(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1],
            );
        }

        /// Matrix3x2-Matrix2x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x2(A: *const Matrix3x2(T), B: *const Matrix2x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1],
                ),
            } };
        }

        /// Matrix3x2-Matrix2x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x3(A: *const Matrix3x2(T), B: *const Matrix2x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1],
                ),
            } };
        }

        /// Matrix3x2-Matrix2x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x4(A: *const Matrix3x2(T), B: *const Matrix2x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1],
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
        pub inline fn div(A: *const Matrix3x2(T), B: *const Matrix3x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).div(&A.m[0], &B.m[0]),
                Vector3(T).div(&A.m[1], &B.m[1]),
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
        pub inline fn divScalar(A: *const Matrix3x2(T), s: T) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).divScalar(&A.m[0], s),
                Vector3(T).divScalar(&A.m[1], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix3x2(T)) Matrix2x3(T) {
            return Matrix2x3(T){ .m = [3]Vector2(T){
                Vector2(T).init(A.m[0].v[0], A.m[1].v[0]),
                Vector2(T).init(A.m[0].v[1], A.m[1].v[1]),
                Vector2(T).init(A.m[0].v[2], A.m[1].v[2]),
            } };
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix3x2(T), B: *const Matrix3x2(T)) bool {
            return Vector3(T).equal(&A.m[0], &B.m[0]) and Vector3(T).equal(&A.m[1], &B.m[1]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix3x2(T), B: *const Matrix3x2(T), tolerance: T) bool {
            return Vector3(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector3(T).approxEqual(&A.m[1], &B.m[1], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        /// - `r2`: The vector for row 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector3(T), r1: *const Vector3(T), r2: *const Vector3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(r0.v[0], r1.v[0], r2.v[0]),
                Vector3(T).init(r0.v[1], r1.v[1], r2.v[1]),
                Vector3(T).init(r0.v[2], r1.v[2], r2.v[2]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        /// - `c2`: The vector for column 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initCol(c0: *const Vector3(T), c1: *const Vector3(T), c2: *const Vector3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(c0.v[0], c0.v[1], c0.v[2]),
                Vector3(T).init(c1.v[0], c1.v[1], c1.v[2]),
                Vector3(T).init(c2.v[0], c2.v[1], c2.v[2]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).splat(s),
                Vector3(T).splat(s),
                Vector3(T).splat(s),
            } };
        }

        /// Identity matrix.
        pub const identity = Matrix3x3(T).init(1, 0, 0, 0, 1, 0, 0, 0, 1);

        /// Zero matrix.
        pub const zero = Matrix3x3(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix3x3(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix3x3(T), r: usize) Vector3(T) {
            return Vector3(T).init(A.m[0].v[r], A.m[1].v[r], A.m[2].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix3x3(T), c: usize) Vector3(T) {
            return Vector3(T).init(A.m[c].v[0], A.m[c].v[1], A.m[c].v[2]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix3x3(T), B: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).add(&A.m[0], &B.m[0]),
                Vector3(T).add(&A.m[1], &B.m[1]),
                Vector3(T).add(&A.m[2], &B.m[2]),
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
        pub inline fn addScalar(A: *const Matrix3x3(T), s: T) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).addScalar(&A.m[0], s),
                Vector3(T).addScalar(&A.m[1], s),
                Vector3(T).addScalar(&A.m[2], s),
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
        pub inline fn sub(A: *const Matrix3x3(T), B: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).sub(&A.m[0], &B.m[0]),
                Vector3(T).sub(&A.m[1], &B.m[1]),
                Vector3(T).sub(&A.m[2], &B.m[2]),
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
        pub inline fn subScalar(A: *const Matrix3x3(T), s: T) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).subScalar(&A.m[0], s),
                Vector3(T).subScalar(&A.m[1], s),
                Vector3(T).subScalar(&A.m[2], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).negate(&A.m[0]),
                Vector3(T).negate(&A.m[1]),
                Vector3(T).negate(&A.m[2]),
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
        pub inline fn mul(A: *const Matrix3x3(T), B: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).mul(&A.m[0], &B.m[0]),
                Vector3(T).mul(&A.m[1], &B.m[1]),
                Vector3(T).mul(&A.m[2], &B.m[2]),
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
        pub inline fn mulScalar(A: *const Matrix3x3(T), s: T) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).mulScalar(&A.m[0], s),
                Vector3(T).mulScalar(&A.m[1], s),
                Vector3(T).mulScalar(&A.m[2], s),
            } };
        }

        /// Matrix3x3-Vector2 multiplication.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector2(A: *const Matrix3x3(T), v: *const Vector2(T)) Vector2(T) {
            return Vector2(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1],
            );
        }

        /// Matrix3x3-Vector3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector3(A: *const Matrix3x3(T), v: *const Vector3(T)) Vector3(T) {
            return Vector3(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1] + A.m[2].v[2] * v.v[2],
            );
        }

        /// Matrix3x3-Matrix3x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x2(A: *const Matrix3x3(T), B: *const Matrix3x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2],
                ),
            } };
        }

        /// Matrix3x3-Matrix3x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x3(A: *const Matrix3x3(T), B: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2],
                ),
            } };
        }

        /// Matrix3x3-Matrix3x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x4(A: *const Matrix3x3(T), B: *const Matrix3x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[2].v[0] * B.m[3].v[2],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[2].v[1] * B.m[3].v[2],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1] + A.m[2].v[2] * B.m[3].v[2],
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
        pub inline fn div(A: *const Matrix3x3(T), B: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).div(&A.m[0], &B.m[0]),
                Vector3(T).div(&A.m[1], &B.m[1]),
                Vector3(T).div(&A.m[2], &B.m[2]),
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
        pub inline fn divScalar(A: *const Matrix3x3(T), s: T) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).divScalar(&A.m[0], s),
                Vector3(T).divScalar(&A.m[1], s),
                Vector3(T).divScalar(&A.m[2], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(A.m[0].v[0], A.m[1].v[0], A.m[2].v[0]),
                Vector3(T).init(A.m[0].v[1], A.m[1].v[1], A.m[2].v[1]),
                Vector3(T).init(A.m[0].v[2], A.m[1].v[2], A.m[2].v[2]),
            } };
        }

        /// Determinant of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The determinant of the matrix.
        pub inline fn determinant(A: *const Matrix3x3(T)) T {
            return A.m[0].v[0] * (A.m[1].v[1] * A.m[2].v[2] - A.m[1].v[2] * A.m[2].v[1]) - A.m[0].v[1] * (A.m[1].v[0] * A.m[2].v[2] - A.m[1].v[2] * A.m[2].v[0]) + A.m[0].v[2] * (A.m[1].v[0] * A.m[2].v[1] - A.m[1].v[1] * A.m[2].v[0]);
        }

        /// Inverse of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The inverse of the matrix.
        pub inline fn inverse(A: *const Matrix3x3(T)) !Matrix3x3(T) {
            const det = Matrix3x3(T).determinant(A);
            if (det == 0) {
                return error.SingularMatrix;
            }

            const invDet = 1 / det;
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(
                    (A.m[1].v[1] * A.m[2].v[2] - A.m[1].v[2] * A.m[2].v[1]) * invDet,
                    (A.m[0].v[2] * A.m[2].v[1] - A.m[0].v[1] * A.m[2].v[2]) * invDet,
                    (A.m[0].v[1] * A.m[1].v[2] - A.m[0].v[2] * A.m[1].v[1]) * invDet,
                ),
                Vector3(T).init(
                    (A.m[1].v[2] * A.m[2].v[0] - A.m[1].v[0] * A.m[2].v[2]) * invDet,
                    (A.m[0].v[0] * A.m[2].v[2] - A.m[0].v[2] * A.m[2].v[0]) * invDet,
                    (A.m[0].v[2] * A.m[1].v[0] - A.m[0].v[0] * A.m[1].v[2]) * invDet,
                ),
                Vector3(T).init(
                    (A.m[1].v[0] * A.m[2].v[1] - A.m[1].v[1] * A.m[2].v[0]) * invDet,
                    (A.m[0].v[1] * A.m[2].v[0] - A.m[0].v[0] * A.m[2].v[1]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[1] - A.m[0].v[1] * A.m[1].v[0]) * invDet,
                ),
            } };
        }

        /// Trace of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The trace of the matrix.
        pub inline fn trace(A: *const Matrix3x3(T)) T {
            return A.m[0].v[0] + A.m[1].v[1] + A.m[2].v[2];
        }

        /// Generate the scaling matrix.
        ///
        /// **Parameters**:
        /// - `x`: The scaling factor along the x-axis.
        /// - `y`: The scaling factor along the y-axis.
        ///
        /// **Returns**:
        /// - The scaling matrix.
        pub inline fn scale(x: T, y: T) Matrix3x3(T) {
            return Matrix3x3(T).init(x, 0, 0, 0, y, 0, 0, 0, 1);
        }

        /// Generate the inverse scaling matrix.
        ///
        /// **Parameters**:
        /// - `A`: The scaling matrix.
        ///
        /// **Returns**:
        /// - The inverse scaling matrix.
        pub inline fn scaleInverse(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T).init(1 / A.m[0].v[0], 0, 0, 0, 1 / A.m[1].v[1], 0, 0, 0, 1);
        }

        /// Generate the shearing matrix along the x-axis.
        ///
        /// **Parameters**:
        /// - `y`: The shearing factor along the y-axis.
        ///
        /// **Returns**:
        /// - The shearing matrix.
        pub inline fn shearX(y: T) Matrix3x3(T) {
            return Matrix3x3(T).init(1, y, 0, 0, 1, 0, 0, 0, 1);
        }

        /// Generate the shearing matrix along the y-axis.
        ///
        /// **Parameters**:
        /// - `x`: The shearing factor along the x-axis.
        ///
        /// **Returns**:
        /// - The shearing matrix.
        pub inline fn shearY(x: T) Matrix3x3(T) {
            return Matrix3x3(T).init(1, 0, 0, x, 1, 0, 0, 0, 1);
        }

        /// Generate the inverse shearing matrix.
        ///
        /// **Parameters**:
        /// - `A`: The shearing matrix.
        ///
        /// **Returns**:
        /// - The inverse shearing matrix.
        pub inline fn shearInverse(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T).init(1, -A.m[1].v[0], 0, -A.m[0].v[1], 1, 0, 0, 0, 1);
        }

        /// Generate the rotation matrix.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotate(angle: T) Matrix3x3(T) {
            const c = @cos(angle);
            const s = @sin(angle);
            return Matrix3x3(T).init(c, s, 0, -s, c, 0, 0, 0, 1);
        }

        /// Generate the inverse rotation matrix.
        ///
        /// **Parameters**:
        /// - `A`: The rotation matrix.
        ///
        /// **Returns**:
        /// - The inverse rotation matrix.
        pub inline fn rotateInverse(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return A.transpose();
        }

        /// Generate the reflection matrix across a normal vector. The
        /// reflection is performed with respect to the line that is
        /// perpendicular to the normal and passes through the origin.
        ///
        /// **Parameters**:
        /// - `n`: The normal vector. Must be normalized.
        ///
        /// **Returns**:
        /// - The reflection matrix.
        pub inline fn reflect(normal: *const Vector2(T)) Matrix3x3(T) {
            //return Matrix2x2(T).init(
            //    1 - 2 * normal.v[0] * normal.v[0],
            //    2 * normal.v[0] * normal.v[1],
            //    2 * normal.v[0] * normal.v[1],
            //    -1 + 2 * normal.v[1] * normal.v[1],
            //);
            return Matrix3x3(T).init(
                1 - 2 * normal.v[0] * normal.v[0],
                -2 * normal.v[0] * normal.v[1],
                0,
                -2 * normal.v[0] * normal.v[1],
                1 - 2 * normal.v[1] * normal.v[1],
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the reflection matrix across a normal vector and a point.
        /// The reflection is performed with respect to the line that passes
        /// through the point and is perpendicular to the normal.
        ///
        /// **Parameters**:
        /// - `n`: The normal vector. Must be normalized.
        /// - `p`: The point.
        ///
        /// **Returns**:
        /// - The reflection matrix.
        pub inline fn reflectPoint(normal: *const Vector2(T), point: *const Vector2(T)) Matrix3x3(T) {
            const d = Vector2(T).dot(normal, point);
            return Matrix3x3(T).init(
                1 - 2 * normal.v[0] * normal.v[0],
                -2 * normal.v[0] * normal.v[1],
                2 * normal.v[0] * d,
                -2 * normal.v[0] * normal.v[1],
                1 - 2 * normal.v[1] * normal.v[1],
                2 * normal.v[1] * d,
                0,
                0,
                1,
            );
        }

        /// Generate the translation matrix.
        ///
        /// **Parameters**:
        /// - `x`: The translation along the x-axis.
        /// - `y`: The translation along the y-axis.
        ///
        /// **Returns**:
        /// - The translation matrix.
        pub inline fn translate(x: T, y: T) Matrix3x3(T) {
            return Matrix3x3(T).init(1, 0, x, 0, 1, y, 0, 0, 1);
        }

        /// Generate the inverse translation matrix.
        ///
        /// **Parameters**:
        /// - `A`: The translation matrix.
        ///
        /// **Returns**:
        /// - The inverse translation matrix.
        pub inline fn translateInverse(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T).init(1, 0, -A.m[2].v[0], 0, 1, -A.m[2].v[1], 0, 0, 1);
        }

        /// Creates the transformation matrix and applies it to the matrix
        /// in-place.
        ///
        /// **Parameters**:
        /// - `A`: The matrix to transform.
        /// - `f`: The transformation matrix generator function.
        /// - `args`: The function arguments.
        ///
        /// **Returns**:
        /// - The transformed matrix.
        pub inline fn applyTransform(A: *Matrix3x3(T), f: anytype, args: anytype) void {
            A.* = @call(.auto, f, args).mulMatrix3x3(A);
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix3x3(T), B: *const Matrix3x3(T)) bool {
            return Vector3(T).equal(&A.m[0], &B.m[0]) and
                Vector3(T).equal(&A.m[1], &B.m[1]) and
                Vector3(T).equal(&A.m[2], &B.m[2]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix3x3(T), B: *const Matrix3x3(T), tolerance: T) bool {
            return Vector3(T).approxEqual(&A.m[0], &B.m[0], tolerance) and
                Vector3(T).approxEqual(&A.m[1], &B.m[1], tolerance) and
                Vector3(T).approxEqual(&A.m[2], &B.m[2], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        /// - `r2`: The vector for row 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector4(T), r1: *const Vector4(T), r2: *const Vector4(T)) Matrix3x4(T) {
            return Matrix3x3(T){ .m = [4]Vector3(T){
                Vector3(T).init(r0.v[0], r1.v[0], r2.v[0]),
                Vector3(T).init(r0.v[1], r1.v[1], r2.v[1]),
                Vector3(T).init(r0.v[2], r1.v[2], r2.v[2]),
                Vector3(T).init(r0.v[3], r1.v[3], r2.v[3]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        /// - `c2`: The vector for column 2.
        /// - `c3`: The vector for column 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initCol(c0: *const Vector3(T), c1: *const Vector3(T), c2: *const Vector3(T), c3: *const Vector3(T)) Matrix3x4(T) {
            return Matrix3x3(T){ .m = [4]Vector3(T){
                Vector3(T).init(c0.v[0], c0.v[1], c0.v[2]),
                Vector3(T).init(c1.v[0], c1.v[1], c1.v[2]),
                Vector3(T).init(c2.v[0], c2.v[1], c2.v[2]),
                Vector3(T).init(c3.v[0], c3.v[1], c3.v[2]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).splat(s),
                Vector3(T).splat(s),
                Vector3(T).splat(s),
                Vector3(T).splat(s),
            } };
        }

        /// Zero matrix.
        pub const zero = Matrix3x4(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix3x4(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix3x4(T), r: usize) Vector4(T) {
            return Vector4(T).init(A.m[0].v[r], A.m[1].v[r], A.m[2].v[r], A.m[3].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix3x4(T), c: usize) Vector3(T) {
            return Vector3(T).init(A.m[c].v[0], A.m[c].v[1], A.m[c].v[2]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix3x4(T), B: *const Matrix3x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).add(&A.m[0], &B.m[0]),
                Vector3(T).add(&A.m[1], &B.m[1]),
                Vector3(T).add(&A.m[2], &B.m[2]),
                Vector3(T).add(&A.m[3], &B.m[3]),
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
        pub inline fn addScalar(A: *const Matrix3x4(T), s: T) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).addScalar(&A.m[0], s),
                Vector3(T).addScalar(&A.m[1], s),
                Vector3(T).addScalar(&A.m[2], s),
                Vector3(T).addScalar(&A.m[3], s),
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
        pub inline fn sub(A: *const Matrix3x4(T), B: *const Matrix3x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).sub(&A.m[0], &B.m[0]),
                Vector3(T).sub(&A.m[1], &B.m[1]),
                Vector3(T).sub(&A.m[2], &B.m[2]),
                Vector3(T).sub(&A.m[3], &B.m[3]),
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
        pub inline fn subScalar(A: *const Matrix3x4(T), s: T) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).subScalar(&A.m[0], s),
                Vector3(T).subScalar(&A.m[1], s),
                Vector3(T).subScalar(&A.m[2], s),
                Vector3(T).subScalar(&A.m[3], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix3x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).negate(&A.m[0]),
                Vector3(T).negate(&A.m[1]),
                Vector3(T).negate(&A.m[2]),
                Vector3(T).negate(&A.m[3]),
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
        pub inline fn mul(A: *const Matrix3x4(T), B: *const Matrix3x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).mul(&A.m[0], &B.m[0]),
                Vector3(T).mul(&A.m[1], &B.m[1]),
                Vector3(T).mul(&A.m[2], &B.m[2]),
                Vector3(T).mul(&A.m[3], &B.m[3]),
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
        pub inline fn mulScalar(A: *const Matrix3x4(T), s: T) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).mulScalar(&A.m[0], s),
                Vector3(T).mulScalar(&A.m[1], s),
                Vector3(T).mulScalar(&A.m[2], s),
                Vector3(T).mulScalar(&A.m[3], s),
            } };
        }

        /// Matrix3x4-Vector4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector4(A: *const Matrix3x4(T), v: *const Vector4(T)) Vector3(T) {
            return Vector3(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2] + A.m[3].v[0] * v.v[3],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2] + A.m[3].v[1] * v.v[3],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1] + A.m[2].v[2] * v.v[2] + A.m[3].v[2] * v.v[3],
            );
        }

        /// Matrix3x4-Matrix4x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x2(A: *const Matrix3x4(T), B: *const Matrix4x2(T)) Matrix3x2(T) {
            return Matrix3x2(T){ .m = [2]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2] + A.m[3].v[2] * B.m[0].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2] + A.m[3].v[2] * B.m[1].v[3],
                ),
            } };
        }

        /// Matrix3x4-Matrix4x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x3(A: *const Matrix3x4(T), B: *const Matrix4x3(T)) Matrix3x3(T) {
            return Matrix3x3(T){ .m = [3]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2] + A.m[3].v[2] * B.m[0].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2] + A.m[3].v[2] * B.m[1].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2] + A.m[3].v[0] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2] + A.m[3].v[1] * B.m[2].v[3],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2] + A.m[3].v[2] * B.m[2].v[3],
                ),
            } };
        }

        /// Matrix3x4-Matrix4x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x4(A: *const Matrix3x4(T), B: *const Matrix4x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2] + A.m[3].v[2] * B.m[0].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2] + A.m[3].v[2] * B.m[1].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2] + A.m[3].v[0] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2] + A.m[3].v[1] * B.m[2].v[3],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2] + A.m[3].v[2] * B.m[2].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[2].v[0] * B.m[3].v[2] + A.m[3].v[0] * B.m[3].v[3],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[2].v[1] * B.m[3].v[2] + A.m[3].v[1] * B.m[3].v[3],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1] + A.m[2].v[2] * B.m[3].v[2] + A.m[3].v[2] * B.m[3].v[3],
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
        pub inline fn div(A: *const Matrix3x4(T), B: *const Matrix3x4(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).div(&A.m[0], &B.m[0]),
                Vector3(T).div(&A.m[1], &B.m[1]),
                Vector3(T).div(&A.m[2], &B.m[2]),
                Vector3(T).div(&A.m[3], &B.m[3]),
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
        pub inline fn divScalar(A: *const Matrix3x4(T), s: T) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).divScalar(&A.m[0], s),
                Vector3(T).divScalar(&A.m[1], s),
                Vector3(T).divScalar(&A.m[2], s),
                Vector3(T).divScalar(&A.m[3], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix3x4(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(A.m[0].v[0], A.m[1].v[0], A.m[2].v[0], A.m[3].v[0]),
                Vector4(T).init(A.m[0].v[1], A.m[1].v[1], A.m[2].v[1], A.m[3].v[1]),
                Vector4(T).init(A.m[0].v[2], A.m[1].v[2], A.m[2].v[2], A.m[3].v[2]),
            } };
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix3x4(T), B: *const Matrix3x4(T)) bool {
            return Vector3(T).equal(&A.m[0], &B.m[0]) and Vector3(T).equal(&A.m[1], &B.m[1]) and Vector3(T).equal(&A.m[2], &B.m[2]) and Vector3(T).equal(&A.m[3], &B.m[3]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix3x4(T), B: *const Matrix3x4(T), tolerance: T) bool {
            return Vector3(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector3(T).approxEqual(&A.m[1], &B.m[1], tolerance) and Vector3(T).approxEqual(&A.m[2], &B.m[2], tolerance) and Vector3(T).approxEqual(&A.m[3], &B.m[3], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        /// - `r2`: The vector for row 2.
        /// - `r3`: The vector for row 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector2(T), r1: *const Vector2(T), r2: *const Vector2(T), r3: *const Vector2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(r0.v[0], r1.v[0], r2.v[0], r3.v[0]),
                Vector4(T).init(r0.v[1], r1.v[1], r2.v[1], r3.v[1]),
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
        pub inline fn initCol(c0: *const Vector4(T), c1: *const Vector4(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(c0.v[0], c0.v[1], c0.v[2], c0.v[3]),
                Vector4(T).init(c1.v[0], c1.v[1], c1.v[2], c1.v[3]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).splat(s),
                Vector4(T).splat(s),
            } };
        }

        /// Zero matrix.
        pub const zero = Matrix4x2(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix4x2(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix4x2(T), r: usize) Vector4(T) {
            return Vector4(T).init(A.m[0].v[r], A.m[1].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix4x2(T), c: usize) Vector2(T) {
            return Vector2(T).init(A.m[c].v[0], A.m[c].v[1]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix4x2(T), B: *const Matrix4x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).add(&A.m[0], &B.m[0]),
                Vector4(T).add(&A.m[1], &B.m[1]),
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
        pub inline fn addScalar(A: *const Matrix4x2(T), s: T) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).addScalar(&A.m[0], s),
                Vector4(T).addScalar(&A.m[1], s),
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
        pub inline fn sub(A: *const Matrix4x2(T), B: *const Matrix4x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).sub(&A.m[0], &B.m[0]),
                Vector4(T).sub(&A.m[1], &B.m[1]),
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
        pub inline fn subScalar(A: *const Matrix4x2(T), s: T) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).subScalar(&A.m[0], s),
                Vector4(T).subScalar(&A.m[1], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix4x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).negate(&A.m[0]),
                Vector4(T).negate(&A.m[1]),
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
        pub inline fn mul(A: *const Matrix4x2(T), B: *const Matrix4x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).mul(&A.m[0], &B.m[0]),
                Vector4(T).mul(&A.m[1], &B.m[1]),
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
        pub inline fn mulScalar(A: *const Matrix4x2(T), s: T) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).mulScalar(&A.m[0], s),
                Vector4(T).mulScalar(&A.m[1], s),
            } };
        }

        /// Matrix4x2-Vector2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector2(A: *const Matrix4x2(T), v: *const Vector2(T)) Vector4(T) {
            return Vector4(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1],
                A.m[0].v[3] * v.v[0] + A.m[1].v[3] * v.v[1],
            );
        }

        /// Matrix4x2-Matrix2x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x2(A: *const Matrix4x2(T), B: *const Matrix2x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1],
                ),
            } };
        }

        /// Matrix4x2-Matrix2x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x3(A: *const Matrix4x2(T), B: *const Matrix2x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1],
                    A.m[0].v[3] * B.m[2].v[0] + A.m[1].v[3] * B.m[2].v[1],
                ),
            } };
        }

        /// Matrix4x2-Matrix2x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix2x4(A: *const Matrix4x2(T), B: *const Matrix2x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1],
                    A.m[0].v[3] * B.m[2].v[0] + A.m[1].v[3] * B.m[2].v[1],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1],
                    A.m[0].v[3] * B.m[3].v[0] + A.m[1].v[3] * B.m[3].v[1],
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
        pub inline fn div(A: *const Matrix4x2(T), B: *const Matrix4x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).div(&A.m[0], &B.m[0]),
                Vector4(T).div(&A.m[1], &B.m[1]),
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
        pub inline fn divScalar(A: *const Matrix4x2(T), s: T) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).divScalar(&A.m[0], s),
                Vector4(T).divScalar(&A.m[1], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix4x2(T)) Matrix2x4(T) {
            return Matrix2x4(T){ .m = [4]Vector2(T){
                Vector2(T).init(A.m[0].v[0], A.m[1].v[0]),
                Vector2(T).init(A.m[0].v[1], A.m[1].v[1]),
                Vector2(T).init(A.m[0].v[2], A.m[1].v[2]),
                Vector2(T).init(A.m[0].v[3], A.m[1].v[3]),
            } };
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix4x2(T), B: *const Matrix4x2(T)) bool {
            return Vector4(T).equal(&A.m[0], &B.m[0]) and Vector4(T).equal(&A.m[1], &B.m[1]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix4x2(T), B: *const Matrix4x2(T), tolerance: T) bool {
            return Vector4(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector4(T).approxEqual(&A.m[1], &B.m[1], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        /// - `r2`: The vector for row 2.
        /// - `r3`: The vector for row 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector3(T), r1: *const Vector3(T), r2: *const Vector3(T), r3: *const Vector3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(r0.v[0], r1.v[0], r2.v[0], r3.v[0]),
                Vector4(T).init(r0.v[1], r1.v[1], r2.v[1], r3.v[1]),
                Vector4(T).init(r0.v[2], r1.v[2], r2.v[2], r3.v[2]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        /// - `c2`: The vector for column 2.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initCol(c0: *const Vector4(T), c1: *const Vector4(T), c2: *const Vector4(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(c0.v[0], c0.v[1], c0.v[2], c0.v[3]),
                Vector4(T).init(c1.v[0], c1.v[1], c1.v[2], c1.v[3]),
                Vector4(T).init(c2.v[0], c2.v[1], c2.v[2], c2.v[3]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).splat(s),
                Vector4(T).splat(s),
                Vector4(T).splat(s),
            } };
        }

        /// Zero matrix.
        pub const zero = Matrix4x3(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix4x3(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix4x3(T), r: usize) Vector4(T) {
            return Vector4(T).init(A.m[0].v[r], A.m[1].v[r], A.m[2].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix4x3(T), c: usize) Vector3(T) {
            return Vector3(T).init(A.m[c].v[0], A.m[c].v[1], A.m[c].v[2]);
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix4x3(T), B: *const Matrix4x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).add(&A.m[0], &B.m[0]),
                Vector4(T).add(&A.m[1], &B.m[1]),
                Vector4(T).add(&A.m[2], &B.m[2]),
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
        pub inline fn addScalar(A: *const Matrix4x3(T), s: T) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).addScalar(&A.m[0], s),
                Vector4(T).addScalar(&A.m[1], s),
                Vector4(T).addScalar(&A.m[2], s),
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
        pub inline fn sub(A: *const Matrix4x3(T), B: *const Matrix4x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).sub(&A.m[0], &B.m[0]),
                Vector4(T).sub(&A.m[1], &B.m[1]),
                Vector4(T).sub(&A.m[2], &B.m[2]),
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
        pub inline fn subScalar(A: *const Matrix4x3(T), s: T) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).subScalar(&A.m[0], s),
                Vector4(T).subScalar(&A.m[1], s),
                Vector4(T).subScalar(&A.m[2], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix4x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).negate(&A.m[0]),
                Vector4(T).negate(&A.m[1]),
                Vector4(T).negate(&A.m[2]),
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
        pub inline fn mul(A: *const Matrix4x3(T), B: *const Matrix4x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).mul(&A.m[0], &B.m[0]),
                Vector4(T).mul(&A.m[1], &B.m[1]),
                Vector4(T).mul(&A.m[2], &B.m[2]),
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
        pub inline fn mulScalar(A: *const Matrix4x3(T), s: T) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).mulScalar(&A.m[0], s),
                Vector4(T).mulScalar(&A.m[1], s),
                Vector4(T).mulScalar(&A.m[2], s),
            } };
        }

        /// Matrix4x3-Vector3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector3(A: *const Matrix4x3(T), v: *const Vector3(T)) Vector4(T) {
            return Vector4(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1] + A.m[2].v[2] * v.v[2],
                A.m[0].v[3] * v.v[0] + A.m[1].v[3] * v.v[1] + A.m[2].v[3] * v.v[2],
            );
        }

        /// Matrix4x3-Matrix3x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x2(A: *const Matrix4x3(T), B: *const Matrix3x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1] + A.m[2].v[3] * B.m[0].v[2],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1] + A.m[2].v[3] * B.m[1].v[2],
                ),
            } };
        }

        /// Matrix4x3-Matrix3x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x3(A: *const Matrix4x3(T), B: *const Matrix3x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1] + A.m[2].v[3] * B.m[0].v[2],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1] + A.m[2].v[3] * B.m[1].v[2],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2],
                    A.m[0].v[3] * B.m[2].v[0] + A.m[1].v[3] * B.m[2].v[1] + A.m[2].v[3] * B.m[2].v[2],
                ),
            } };
        }

        /// Matrix4x3-Matrix3x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix3x4(A: *const Matrix4x3(T), B: *const Matrix3x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1] + A.m[2].v[3] * B.m[0].v[2],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1] + A.m[2].v[3] * B.m[1].v[2],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2],
                    A.m[0].v[3] * B.m[2].v[0] + A.m[1].v[3] * B.m[2].v[1] + A.m[2].v[3] * B.m[2].v[2],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[2].v[0] * B.m[3].v[2],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[2].v[1] * B.m[3].v[2],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1] + A.m[2].v[2] * B.m[3].v[2],
                    A.m[0].v[3] * B.m[3].v[0] + A.m[1].v[3] * B.m[3].v[1] + A.m[2].v[3] * B.m[3].v[2],
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
        pub inline fn div(A: *const Matrix4x3(T), B: *const Matrix4x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).div(&A.m[0], &B.m[0]),
                Vector4(T).div(&A.m[1], &B.m[1]),
                Vector4(T).div(&A.m[2], &B.m[2]),
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
        pub inline fn divScalar(A: *const Matrix4x3(T), s: T) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).divScalar(&A.m[0], s),
                Vector4(T).divScalar(&A.m[1], s),
                Vector4(T).divScalar(&A.m[2], s),
            } };
        }

        /// Transpose a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix4x3(T)) Matrix3x4(T) {
            return Matrix3x4(T){ .m = [4]Vector3(T){
                Vector3(T).init(A.m[0].v[0], A.m[1].v[0], A.m[2].v[0]),
                Vector3(T).init(A.m[0].v[1], A.m[1].v[1], A.m[2].v[1]),
                Vector3(T).init(A.m[0].v[2], A.m[1].v[2], A.m[2].v[2]),
                Vector3(T).init(A.m[0].v[3], A.m[1].v[3], A.m[2].v[3]),
            } };
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix4x3(T), B: *const Matrix4x3(T)) bool {
            return Vector4(T).equal(&A.m[0], &B.m[0]) and Vector4(T).equal(&A.m[1], &B.m[1]) and Vector4(T).equal(&A.m[2], &B.m[2]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix4x3(T), B: *const Matrix4x3(T), tolerance: T) bool {
            return Vector4(T).approxEqual(&A.m[0], &B.m[0], tolerance) and Vector4(T).approxEqual(&A.m[1], &B.m[1], tolerance) and Vector4(T).approxEqual(&A.m[2], &B.m[2], tolerance);
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

        /// Initialize a matrix with the given row vectors.
        ///
        /// **Parameters**:
        /// - `r0`: The vector for row 0.
        /// - `r1`: The vector for row 1.
        /// - `r2`: The vector for row 2.
        /// - `r3`: The vector for row 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given row vectors.
        pub inline fn initRow(r0: *const Vector4(T), r1: *const Vector4(T), r2: *const Vector4(T), r3: *const Vector4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(r0.v[0], r1.v[0], r2.v[0], r3.v[0]),
                Vector4(T).init(r0.v[1], r1.v[1], r2.v[1], r3.v[1]),
                Vector4(T).init(r0.v[2], r1.v[2], r2.v[2], r3.v[2]),
                Vector4(T).init(r0.v[3], r1.v[3], r2.v[3], r3.v[3]),
            } };
        }

        /// Initialize a matrix with the given column vectors.
        ///
        /// **Parameters**:
        /// - `c0`: The vector for column 0.
        /// - `c1`: The vector for column 1.
        /// - `c2`: The vector for column 2.
        /// - `c3`: The vector for column 3.
        ///
        /// **Returns**:
        /// - A new matrix with the given column vectors.
        pub inline fn initCol(c0: *const Vector4(T), c1: *const Vector4(T), c2: *const Vector4(T), c3: *const Vector4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(c0.v[0], c0.v[1], c0.v[2], c0.v[3]),
                Vector4(T).init(c1.v[0], c1.v[1], c1.v[2], c1.v[3]),
                Vector4(T).init(c2.v[0], c2.v[1], c2.v[2], c2.v[3]),
                Vector4(T).init(c3.v[0], c3.v[1], c3.v[2], c3.v[3]),
            } };
        }

        /// Creates a new matrix with the given scalar value in all components.
        ///
        /// **Parameters**:
        /// - `s`: The scalar value to splat.
        ///
        /// **Returns**:
        /// - A new matrix with the given scalar value in all components.
        pub inline fn splat(s: T) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).splat(s),
                Vector4(T).splat(s),
                Vector4(T).splat(s),
                Vector4(T).splat(s),
            } };
        }

        /// Identity matrix.
        pub const identity = Matrix4x4(T).init(
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            1,
        );

        /// Zero matrix.
        pub const zero = Matrix4x4(T).splat(0);

        /// Get the value at the given row and column.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `r`: The row index.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The value at the given row and column.
        pub inline fn get(A: *const Matrix4x4(T), r: usize, c: usize) T {
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
        pub inline fn row(A: *const Matrix4x4(T), r: usize) Vector4(T) {
            return Vector4(T).init(A.m[0].v[r], A.m[1].v[r], A.m[2].v[r], A.m[3].v[r]);
        }

        /// Get the column vector at the given index.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `c`: The column index.
        ///
        /// **Returns**:
        /// - The column vector at the given index.
        pub inline fn col(A: *const Matrix4x4(T), c: usize) Vector4(T) {
            return A.m[c];
        }

        /// Element-wise addition.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the element-wise sum of the input matrices.
        pub inline fn add(A: *const Matrix4x4(T), B: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).add(&A.m[0], &B.m[0]),
                Vector4(T).add(&A.m[1], &B.m[1]),
                Vector4(T).add(&A.m[2], &B.m[2]),
                Vector4(T).add(&A.m[3], &B.m[3]),
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
        pub inline fn addScalar(A: *const Matrix4x4(T), s: T) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).addScalar(&A.m[0], s),
                Vector4(T).addScalar(&A.m[1], s),
                Vector4(T).addScalar(&A.m[2], s),
                Vector4(T).addScalar(&A.m[3], s),
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
        pub inline fn sub(A: *const Matrix4x4(T), B: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).sub(&A.m[0], &B.m[0]),
                Vector4(T).sub(&A.m[1], &B.m[1]),
                Vector4(T).sub(&A.m[2], &B.m[2]),
                Vector4(T).sub(&A.m[3], &B.m[3]),
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
        pub inline fn subScalar(A: *const Matrix4x4(T), s: T) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).subScalar(&A.m[0], s),
                Vector4(T).subScalar(&A.m[1], s),
                Vector4(T).subScalar(&A.m[2], s),
                Vector4(T).subScalar(&A.m[3], s),
            } };
        }

        /// Negates all components of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - A new matrix with the negated components.
        pub inline fn negate(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).negate(&A.m[0]),
                Vector4(T).negate(&A.m[1]),
                Vector4(T).negate(&A.m[2]),
                Vector4(T).negate(&A.m[3]),
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
        pub inline fn mul(A: *const Matrix4x4(T), B: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).mul(&A.m[0], &B.m[0]),
                Vector4(T).mul(&A.m[1], &B.m[1]),
                Vector4(T).mul(&A.m[2], &B.m[2]),
                Vector4(T).mul(&A.m[3], &B.m[3]),
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
        pub inline fn mulScalar(A: *const Matrix4x4(T), s: T) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).mulScalar(&A.m[0], s),
                Vector4(T).mulScalar(&A.m[1], s),
                Vector4(T).mulScalar(&A.m[2], s),
                Vector4(T).mulScalar(&A.m[3], s),
            } };
        }

        /// Matrix4x4-Vector3 multiplication.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector3(A: *const Matrix4x4(T), v: *const Vector3(T)) Vector3(T) {
            return Vector3(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2] + A.m[3].v[0],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2] + A.m[3].v[1],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1] + A.m[2].v[2] * v.v[2] + A.m[3].v[2],
            );
        }

        /// Matrix4x4-Vector4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{v}`
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        /// - `v`: The vector.
        ///
        /// **Returns**:
        /// - The product of the matrix and vector.
        pub inline fn mulVector4(A: *const Matrix4x4(T), v: *const Vector4(T)) Vector4(T) {
            return Vector4(T).init(
                A.m[0].v[0] * v.v[0] + A.m[1].v[0] * v.v[1] + A.m[2].v[0] * v.v[2] + A.m[3].v[0] * v.v[3],
                A.m[0].v[1] * v.v[0] + A.m[1].v[1] * v.v[1] + A.m[2].v[1] * v.v[2] + A.m[3].v[1] * v.v[3],
                A.m[0].v[2] * v.v[0] + A.m[1].v[2] * v.v[1] + A.m[2].v[2] * v.v[2] + A.m[3].v[2] * v.v[3],
                A.m[0].v[3] * v.v[0] + A.m[1].v[3] * v.v[1] + A.m[2].v[3] * v.v[2] + A.m[3].v[3] * v.v[3],
            );
        }

        /// Matrix4x4-Matrix4x2 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x2(A: *const Matrix4x4(T), B: *const Matrix4x2(T)) Matrix4x2(T) {
            return Matrix4x2(T){ .m = [2]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2] + A.m[3].v[2] * B.m[0].v[3],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1] + A.m[2].v[3] * B.m[0].v[2] + A.m[3].v[3] * B.m[0].v[3],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2] + A.m[3].v[2] * B.m[1].v[3],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1] + A.m[2].v[3] * B.m[1].v[2] + A.m[3].v[3] * B.m[1].v[3],
                ),
            } };
        }

        /// Matrix4x4-Matrix4x3 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x3(A: *const Matrix4x4(T), B: *const Matrix4x3(T)) Matrix4x3(T) {
            return Matrix4x3(T){ .m = [3]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2] + A.m[3].v[2] * B.m[0].v[3],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1] + A.m[2].v[3] * B.m[0].v[2] + A.m[3].v[3] * B.m[0].v[3],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2] + A.m[3].v[2] * B.m[1].v[3],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1] + A.m[2].v[3] * B.m[1].v[2] + A.m[3].v[3] * B.m[1].v[3],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2] + A.m[3].v[0] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2] + A.m[3].v[1] * B.m[2].v[3],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2] + A.m[3].v[2] * B.m[2].v[3],
                    A.m[0].v[3] * B.m[2].v[0] + A.m[1].v[3] * B.m[2].v[1] + A.m[2].v[3] * B.m[2].v[2] + A.m[3].v[3] * B.m[2].v[3],
                ),
            } };
        }

        /// Matrix4x4-Matrix4x4 multiplication:
        ///
        /// :math:`\mathbf{A} \mathbf{B}`
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - The product of the two matrices.
        pub inline fn mulMatrix4x4(A: *const Matrix4x4(T), B: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[2].v[0] * B.m[0].v[2] + A.m[3].v[0] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[2].v[1] * B.m[0].v[2] + A.m[3].v[1] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[2].v[2] * B.m[0].v[2] + A.m[3].v[2] * B.m[0].v[3],
                    A.m[0].v[3] * B.m[0].v[0] + A.m[1].v[3] * B.m[0].v[1] + A.m[2].v[3] * B.m[0].v[2] + A.m[3].v[3] * B.m[0].v[3],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[2].v[0] * B.m[1].v[2] + A.m[3].v[0] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[2].v[1] * B.m[1].v[2] + A.m[3].v[1] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[2].v[2] * B.m[1].v[2] + A.m[3].v[2] * B.m[1].v[3],
                    A.m[0].v[3] * B.m[1].v[0] + A.m[1].v[3] * B.m[1].v[1] + A.m[2].v[3] * B.m[1].v[2] + A.m[3].v[3] * B.m[1].v[3],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[2].v[0] * B.m[2].v[2] + A.m[3].v[0] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[2].v[1] * B.m[2].v[2] + A.m[3].v[1] * B.m[2].v[3],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[2].v[2] * B.m[2].v[2] + A.m[3].v[2] * B.m[2].v[3],
                    A.m[0].v[3] * B.m[2].v[0] + A.m[1].v[3] * B.m[2].v[1] + A.m[2].v[3] * B.m[2].v[2] + A.m[3].v[3] * B.m[2].v[3],
                ),
                Vector4(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[2].v[0] * B.m[3].v[2] + A.m[3].v[0] * B.m[3].v[3],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[2].v[1] * B.m[3].v[2] + A.m[3].v[1] * B.m[3].v[3],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1] + A.m[2].v[2] * B.m[3].v[2] + A.m[3].v[2] * B.m[3].v[3],
                    A.m[0].v[3] * B.m[3].v[0] + A.m[1].v[3] * B.m[3].v[1] + A.m[2].v[3] * B.m[3].v[2] + A.m[3].v[3] * B.m[3].v[3],
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
        pub inline fn div(A: *const Matrix4x4(T), B: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).div(&A.m[0], &B.m[0]),
                Vector4(T).div(&A.m[1], &B.m[1]),
                Vector4(T).div(&A.m[2], &B.m[2]),
                Vector4(T).div(&A.m[3], &B.m[3]),
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
        pub inline fn divScalar(A: *const Matrix4x4(T), s: T) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).divScalar(&A.m[0], s),
                Vector4(T).divScalar(&A.m[1], s),
                Vector4(T).divScalar(&A.m[2], s),
                Vector4(T).divScalar(&A.m[3], s),
            } };
        }

        /// Transposes a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The transposed matrix.
        pub inline fn transpose(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(A.m[0].v[0], A.m[1].v[0], A.m[2].v[0], A.m[3].v[0]),
                Vector4(T).init(A.m[0].v[1], A.m[1].v[1], A.m[2].v[1], A.m[3].v[1]),
                Vector4(T).init(A.m[0].v[2], A.m[1].v[2], A.m[2].v[2], A.m[3].v[2]),
                Vector4(T).init(A.m[0].v[3], A.m[1].v[3], A.m[2].v[3], A.m[3].v[3]),
            } };
        }

        /// Determinant of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The determinant of the matrix.
        pub inline fn determinant(A: *const Matrix4x4(T)) T {
            return A.m[0].v[0] * A.m[1].v[1] * A.m[2].v[2] * A.m[3].v[3] - A.m[0].v[0] * A.m[1].v[1] * A.m[2].v[3] * A.m[3].v[2] +
                A.m[0].v[0] * A.m[1].v[2] * A.m[2].v[3] * A.m[3].v[1] - A.m[0].v[0] * A.m[1].v[2] * A.m[2].v[1] * A.m[3].v[3] +
                A.m[0].v[0] * A.m[1].v[3] * A.m[2].v[1] * A.m[3].v[2] - A.m[0].v[0] * A.m[1].v[3] * A.m[2].v[2] * A.m[3].v[1] -
                A.m[0].v[1] * A.m[1].v[2] * A.m[2].v[3] * A.m[3].v[0] + A.m[0].v[1] * A.m[1].v[2] * A.m[2].v[0] * A.m[3].v[3] -
                A.m[0].v[1] * A.m[1].v[3] * A.m[2].v[0] * A.m[3].v[2] + A.m[0].v[1] * A.m[1].v[3] * A.m[2].v[2] * A.m[3].v[0] -
                A.m[0].v[1] * A.m[1].v[0] * A.m[2].v[2] * A.m[3].v[3] + A.m[0].v[1] * A.m[1].v[0] * A.m[2].v[3] * A.m[3].v[2] +
                A.m[0].v[2] * A.m[1].v[3] * A.m[2].v[0] * A.m[3].v[1] - A.m[0].v[2] * A.m[1].v[3] * A.m[2].v[1] * A.m[3].v[0] +
                A.m[0].v[2] * A.m[1].v[0] * A.m[2].v[1] * A.m[3].v[3] - A.m[0].v[2] * A.m[1].v[0] * A.m[2].v[3] * A.m[3].v[1] +
                A.m[0].v[2] * A.m[1].v[1] * A.m[2].v[3] * A.m[3].v[0] - A.m[0].v[2] * A.m[1].v[1] * A.m[2].v[0] * A.m[3].v[3] -
                A.m[0].v[3] * A.m[1].v[0] * A.m[2].v[1] * A.m[3].v[2] + A.m[0].v[3] * A.m[1].v[0] * A.m[2].v[2] * A.m[3].v[1] -
                A.m[0].v[3] * A.m[1].v[1] * A.m[2].v[2] * A.m[3].v[0] + A.m[0].v[3] * A.m[1].v[1] * A.m[2].v[0] * A.m[3].v[2] -
                A.m[0].v[3] * A.m[1].v[2] * A.m[2].v[0] * A.m[3].v[1] + A.m[0].v[3] * A.m[1].v[2] * A.m[2].v[1] * A.m[3].v[0];
        }

        /// Inverse of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The inverse of the matrix.
        pub inline fn inverse(A: *const Matrix4x4(T)) !Matrix4x4(T) {
            const det = Matrix4x4(T).determinant(A);
            if (det == 0) {
                return error.SingularMatrix;
            }

            const invDet = 1 / det;
            return Matrix4x4(T){ .m = [4]Vector4(T){
                Vector4(T).init(
                    (A.m[1].v[1] * A.m[2].v[2] * A.m[3].v[3] + A.m[1].v[2] * A.m[2].v[3] * A.m[3].v[1] + A.m[1].v[3] * A.m[2].v[1] * A.m[3].v[2] - A.m[1].v[1] * A.m[2].v[3] * A.m[3].v[2] - A.m[1].v[2] * A.m[2].v[1] * A.m[3].v[3] - A.m[1].v[3] * A.m[2].v[2] * A.m[3].v[1]) * invDet,
                    (A.m[0].v[1] * A.m[2].v[3] * A.m[3].v[2] + A.m[0].v[2] * A.m[2].v[1] * A.m[3].v[3] + A.m[0].v[3] * A.m[2].v[2] * A.m[3].v[1] - A.m[0].v[1] * A.m[2].v[2] * A.m[3].v[3] - A.m[0].v[2] * A.m[2].v[3] * A.m[3].v[1] - A.m[0].v[3] * A.m[2].v[1] * A.m[3].v[2]) * invDet,
                    (A.m[0].v[1] * A.m[1].v[2] * A.m[3].v[3] + A.m[0].v[2] * A.m[1].v[3] * A.m[3].v[1] + A.m[0].v[3] * A.m[1].v[1] * A.m[3].v[2] - A.m[0].v[1] * A.m[1].v[3] * A.m[3].v[2] - A.m[0].v[2] * A.m[1].v[1] * A.m[3].v[3] - A.m[0].v[3] * A.m[1].v[2] * A.m[3].v[1]) * invDet,
                    (A.m[0].v[1] * A.m[1].v[3] * A.m[2].v[2] + A.m[0].v[2] * A.m[1].v[1] * A.m[2].v[3] + A.m[0].v[3] * A.m[1].v[2] * A.m[2].v[1] - A.m[0].v[1] * A.m[1].v[2] * A.m[2].v[3] - A.m[0].v[2] * A.m[1].v[3] * A.m[2].v[1] - A.m[0].v[3] * A.m[1].v[1] * A.m[2].v[2]) * invDet,
                ),
                Vector4(T).init(
                    (A.m[1].v[0] * A.m[2].v[3] * A.m[3].v[2] + A.m[1].v[2] * A.m[2].v[0] * A.m[3].v[3] + A.m[1].v[3] * A.m[2].v[2] * A.m[3].v[0] - A.m[1].v[0] * A.m[2].v[2] * A.m[3].v[3] - A.m[1].v[2] * A.m[2].v[3] * A.m[3].v[0] - A.m[1].v[3] * A.m[2].v[0] * A.m[3].v[2]) * invDet,
                    (A.m[0].v[0] * A.m[2].v[2] * A.m[3].v[3] + A.m[0].v[2] * A.m[2].v[3] * A.m[3].v[0] + A.m[0].v[3] * A.m[2].v[0] * A.m[3].v[2] - A.m[0].v[0] * A.m[2].v[3] * A.m[3].v[2] - A.m[0].v[2] * A.m[2].v[0] * A.m[3].v[3] - A.m[0].v[3] * A.m[2].v[2] * A.m[3].v[0]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[3] * A.m[3].v[2] + A.m[0].v[2] * A.m[1].v[0] * A.m[3].v[3] + A.m[0].v[3] * A.m[1].v[2] * A.m[3].v[0] - A.m[0].v[0] * A.m[1].v[2] * A.m[3].v[3] - A.m[0].v[2] * A.m[1].v[3] * A.m[3].v[0] - A.m[0].v[3] * A.m[1].v[0] * A.m[3].v[2]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[2] * A.m[2].v[3] + A.m[0].v[2] * A.m[1].v[3] * A.m[2].v[0] + A.m[0].v[3] * A.m[1].v[0] * A.m[2].v[2] - A.m[0].v[0] * A.m[1].v[3] * A.m[2].v[2] - A.m[0].v[2] * A.m[1].v[0] * A.m[2].v[3] - A.m[0].v[3] * A.m[1].v[2] * A.m[2].v[0]) * invDet,
                ),
                Vector4(T).init(
                    (A.m[1].v[0] * A.m[2].v[1] * A.m[3].v[3] + A.m[1].v[1] * A.m[2].v[3] * A.m[3].v[0] + A.m[1].v[3] * A.m[2].v[0] * A.m[3].v[1] - A.m[1].v[0] * A.m[2].v[3] * A.m[3].v[1] - A.m[1].v[1] * A.m[2].v[0] * A.m[3].v[3] - A.m[1].v[3] * A.m[2].v[1] * A.m[3].v[0]) * invDet,
                    (A.m[0].v[0] * A.m[2].v[3] * A.m[3].v[1] + A.m[0].v[1] * A.m[2].v[0] * A.m[3].v[3] + A.m[0].v[3] * A.m[2].v[1] * A.m[3].v[0] - A.m[0].v[0] * A.m[2].v[1] * A.m[3].v[3] - A.m[0].v[1] * A.m[2].v[3] * A.m[3].v[0] - A.m[0].v[3] * A.m[2].v[0] * A.m[3].v[1]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[1] * A.m[3].v[3] + A.m[0].v[1] * A.m[1].v[3] * A.m[3].v[0] + A.m[0].v[3] * A.m[1].v[0] * A.m[3].v[1] - A.m[0].v[0] * A.m[1].v[3] * A.m[3].v[1] - A.m[0].v[1] * A.m[1].v[0] * A.m[3].v[3] - A.m[0].v[3] * A.m[1].v[1] * A.m[3].v[0]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[3] * A.m[2].v[1] + A.m[0].v[1] * A.m[1].v[0] * A.m[2].v[3] + A.m[0].v[3] * A.m[1].v[1] * A.m[2].v[0] - A.m[0].v[0] * A.m[1].v[1] * A.m[2].v[3] - A.m[0].v[1] * A.m[1].v[3] * A.m[2].v[0] - A.m[0].v[3] * A.m[1].v[0] * A.m[2].v[1]) * invDet,
                ),
                Vector4(T).init(
                    (A.m[1].v[0] * A.m[2].v[2] * A.m[3].v[1] + A.m[1].v[1] * A.m[2].v[0] * A.m[3].v[2] + A.m[1].v[2] * A.m[2].v[1] * A.m[3].v[0] - A.m[1].v[0] * A.m[2].v[1] * A.m[3].v[2] - A.m[1].v[1] * A.m[2].v[2] * A.m[3].v[0] - A.m[1].v[2] * A.m[2].v[0] * A.m[3].v[1]) * invDet,
                    (A.m[0].v[0] * A.m[2].v[1] * A.m[3].v[2] + A.m[0].v[1] * A.m[2].v[2] * A.m[3].v[0] + A.m[0].v[2] * A.m[2].v[0] * A.m[3].v[1] - A.m[0].v[0] * A.m[2].v[2] * A.m[3].v[1] - A.m[0].v[1] * A.m[2].v[0] * A.m[3].v[2] - A.m[0].v[2] * A.m[2].v[1] * A.m[3].v[0]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[2] * A.m[3].v[1] + A.m[0].v[1] * A.m[1].v[0] * A.m[3].v[2] + A.m[0].v[2] * A.m[1].v[1] * A.m[3].v[0] - A.m[0].v[0] * A.m[1].v[1] * A.m[3].v[2] - A.m[0].v[1] * A.m[1].v[2] * A.m[3].v[0] - A.m[0].v[2] * A.m[1].v[0] * A.m[3].v[1]) * invDet,
                    (A.m[0].v[0] * A.m[1].v[1] * A.m[2].v[2] + A.m[0].v[1] * A.m[1].v[2] * A.m[2].v[0] + A.m[0].v[2] * A.m[1].v[0] * A.m[2].v[1] - A.m[0].v[0] * A.m[1].v[2] * A.m[2].v[1] - A.m[0].v[1] * A.m[1].v[0] * A.m[2].v[2] - A.m[0].v[2] * A.m[1].v[1] * A.m[2].v[0]) * invDet,
                ),
            } };
        }

        /// Trace of a matrix.
        ///
        /// **Parameters**:
        /// - `A`: The matrix.
        ///
        /// **Returns**:
        /// - The trace of the matrix.
        pub inline fn trace(A: *const Matrix4x4(T)) T {
            return A.m[0].v[0] + A.m[1].v[1] + A.m[2].v[2] + A.m[3].v[3];
        }

        /// Generate the scaling matrix.
        ///
        /// **Parameters**:
        /// - `x`: The scaling factor along the x-axis.
        /// - `y`: The scaling factor along the y-axis.
        /// - `z`: The scaling factor along the z-axis.
        ///
        /// **Returns**:
        /// - The scaling matrix.
        pub inline fn scale(x: T, y: T, z: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                x,
                0,
                0,
                0,
                0,
                y,
                0,
                0,
                0,
                0,
                z,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the inverse scaling matrix.
        ///
        /// **Parameters**:
        /// - `A`: The scaling matrix.
        ///
        /// **Returns**:
        /// - The inverse scaling matrix.
        pub inline fn scaleInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1 / A.m[0].v[0],
                0,
                0,
                0,
                0,
                1 / A.m[1].v[1],
                0,
                0,
                0,
                0,
                1 / A.m[2].v[2],
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the shearing matrix along the x-axis.
        ///
        /// **Parameters**:
        /// - `y`: The shearing factor along the y-axis.
        /// - `z`: The shearing factor along the z-axis.
        ///
        /// **Returns**:
        /// - The shearing matrix.
        pub inline fn shearX(y: T, z: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1,
                y,
                z,
                0,
                0,
                1,
                0,
                0,
                0,
                0,
                1,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the shearing matrix along the y-axis.
        ///
        /// **Parameters**:
        /// - `x`: The shearing factor along the x-axis.
        /// - `z`: The shearing factor along the z-axis.
        ///
        /// **Returns**:
        /// - The shearing matrix.
        pub inline fn shearY(x: T, z: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1,
                0,
                0,
                0,
                x,
                1,
                z,
                0,
                0,
                0,
                1,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the shearing matrix along the z-axis.
        ///
        /// **Parameters**:
        /// - `x`: The shearing factor along the x-axis.
        /// - `y`: The shearing factor along the y-axis.
        ///
        /// **Returns**:
        /// - The shearing matrix.
        pub inline fn shearZ(x: T, y: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1,
                0,
                0,
                0,
                0,
                1,
                0,
                0,
                x,
                y,
                1,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the inverse shearing matrix.
        ///
        /// **Parameters**:
        /// - `A`: The shearing matrix.
        ///
        /// **Returns**:
        /// - The inverse shearing matrix.
        pub inline fn shearInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1,
                -A.m[1].v[0] / A.m[1].v[1],
                -A.m[2].v[0] / A.m[2].v[2],
                0,
                -A.m[0].v[1] / A.m[0].v[0],
                1,
                -A.m[2].v[1] / A.m[2].v[2],
                0,
                -A.m[0].v[2] / A.m[0].v[0],
                -A.m[1].v[2] / A.m[1].v[1],
                1,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the rotation matrix around an arbitrary axis.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        /// - `axis`: The rotation axis. Must be normalized.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotate(angle: T, axis: *const Vector3(T)) Matrix4x4(T) {
            const c = @cos(angle);
            const s = @sin(angle);
            const t = 1 - c;
            const x = axis.v[0];
            const y = axis.v[1];
            const z = axis.v[2];
            return Matrix4x4(T).init(
                t * x * x + c,
                t * x * y - s * z,
                t * x * z + s * y,
                0,
                t * x * y + s * z,
                t * y * y + c,
                t * y * z - s * x,
                0,
                t * x * z - s * y,
                t * y * z + s * x,
                t * z * z + c,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the rotation matrix around the x-axis.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotateX(angle: T) Matrix4x4(T) {
            const c = @cos(angle);
            const s = @sin(angle);
            return Matrix4x4(T).init(
                1,
                0,
                0,
                0,
                0,
                c,
                -s,
                0,
                0,
                s,
                c,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the rotation matrix around the y-axis.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotateY(angle: T) Matrix4x4(T) {
            const c = @cos(angle);
            const s = @sin(angle);
            return Matrix4x4(T).init(
                c,
                0,
                s,
                0,
                0,
                1,
                0,
                0,
                -s,
                0,
                c,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the rotation matrix around the z-axis.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotateZ(angle: T) Matrix4x4(T) {
            const c = @cos(angle);
            const s = @sin(angle);
            return Matrix4x4(T).init(
                c,
                -s,
                0,
                0,
                s,
                c,
                0,
                0,
                0,
                0,
                1,
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the inverse rotation matrix.
        ///
        /// **Parameters**:
        /// - `A`: The rotation matrix.
        ///
        /// **Returns**:
        /// - The inverse rotation matrix.
        pub inline fn rotateInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).transpose(A);
        }

        /// Generate the reflection matrix across a normal vector. The
        /// reflection is performed with respect to the plane that is
        /// perpendicular to the normal and passes through the origin.
        ///
        /// **Parameters**:
        /// - `n`: The normal vector. Must be normalized.
        ///
        /// **Returns**:
        /// - The reflection matrix.
        pub inline fn reflect(n: *const Vector3(T)) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1 - 2 * n.v[0] * n.v[0],
                -2 * n.v[0] * n.v[1],
                -2 * n.v[0] * n.v[2],
                0,
                -2 * n.v[1] * n.v[0],
                1 - 2 * n.v[1] * n.v[1],
                -2 * n.v[1] * n.v[2],
                0,
                -2 * n.v[2] * n.v[0],
                -2 * n.v[2] * n.v[1],
                1 - 2 * n.v[2] * n.v[2],
                0,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the reflection matrix across a normal vector and a point.
        /// The reflection is performed with respect to the plane that passes
        /// through the point and is perpendicular to the normal.
        ///
        /// **Parameters**:
        /// - `n`: The normal vector. Must be normalized.
        /// - `p`: The point.
        ///
        /// **Returns**:
        /// - The reflection matrix.
        pub inline fn reflectPoint(n: *const Vector3(T), p: *const Vector3(T)) Matrix4x4(T) {
            const d = Vector3(T).dot(n, p);
            return Matrix4x4(T).init(
                1 - 2 * n.v[0] * n.v[0],
                -2 * n.v[0] * n.v[1],
                -2 * n.v[0] * n.v[2],
                2 * d * n.v[0],
                -2 * n.v[1] * n.v[0],
                1 - 2 * n.v[1] * n.v[1],
                -2 * n.v[1] * n.v[2],
                2 * d * n.v[1],
                -2 * n.v[2] * n.v[0],
                -2 * n.v[2] * n.v[1],
                1 - 2 * n.v[2] * n.v[2],
                2 * d * n.v[2],
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the translation matrix.
        ///
        /// **Parameters**:
        /// - `x`: The translation along the x-axis.
        /// - `y`: The translation along the y-axis.
        /// - `z`: The translation along the z-axis.
        ///
        /// **Returns**:
        /// - The translation matrix.
        pub inline fn translate(x: T, y: T, z: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1,
                0,
                0,
                x,
                0,
                1,
                0,
                y,
                0,
                0,
                1,
                z,
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the inverse translation matrix.
        ///
        /// **Parameters**:
        /// - `A`: The translation matrix.
        ///
        /// **Returns**:
        /// - The inverse translation matrix.
        pub inline fn translateInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1,
                0,
                0,
                -A.m[3].v[0],
                0,
                1,
                0,
                -A.m[3].v[1],
                0,
                0,
                1,
                -A.m[3].v[2],
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the view matrix for a coordinate system. The view matrix
        /// transforms from world space to camera space.
        ///
        /// **Parameters**:
        /// - `eye`: The position of the camera.
        /// - `center`: The position to look at.
        /// - `up`: The up direction.
        ///
        /// **Returns**:
        /// - The view matrix.
        pub inline fn lookAt(eye: *const Vector3(T), center: *const Vector3(T), up: *const Vector3(T)) Matrix4x4(T) {
            const f = Vector3(T).normalize(&Vector3(T).sub(center, eye));
            const s = Vector3(T).normalize(&Vector3(T).cross(&f, up));
            const u = Vector3(T).cross(&s, &f);
            return Matrix4x4(T).init(
                s.v[0],
                s.v[1],
                s.v[2],
                -Vector3(T).dot(&s, eye),
                u.v[0],
                u.v[1],
                u.v[2],
                -Vector3(T).dot(&u, eye),
                -f.v[0],
                -f.v[1],
                -f.v[2],
                Vector3(T).dot(&f, eye),
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the inverse view matrix.
        ///
        /// **Parameters**:
        /// - `A`: The view matrix.
        ///
        /// **Returns**:
        /// - The inverse view matrix.
        pub inline fn lookAtInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).transpose(A);
        }

        /// Generate the perspective projection matrix for a negative-one-to-one
        /// depth range. The projection matrix transforms from camera space to
        /// clip space.
        ///
        /// **Parameters**:
        /// - `fov`: The field of view in radians.
        /// - `aspect`: The aspect ratio.
        /// - `near`: The distance to the near clipping plane.
        /// - `far`: The distance to the far clipping plane.
        ///
        /// **Returns**:
        /// - The perspective projection matrix.
        pub inline fn perspectiveNO(fov: T, aspect: T, near: T, far: T) Matrix4x4(T) {
            const f = 1 / @tan(fov / 2);
            return Matrix4x4(T).init(
                f / aspect,
                0,
                0,
                0,
                0,
                f,
                0,
                0,
                0,
                0,
                (far + near) / (near - far),
                (2 * near * far) / (near - far),
                0,
                0,
                -1,
                0,
            );
        }

        /// Generate the perspective projection matrix for a zero-to-one depth
        /// range. The projection matrix transforms from camera space to clip
        /// space.
        ///
        /// **Parameters**:
        /// - `fov`: The field of view in radians.
        /// - `aspect`: The aspect ratio.
        /// - `near`: The distance to the near clipping plane.
        /// - `far`: The distance to the far clipping plane.
        ///
        /// **Returns**:
        /// - The perspective projection matrix.
        pub inline fn perspectiveZO(fov: T, aspect: T, near: T, far: T) Matrix4x4(T) {
            const f = 1 / @tan(fov / 2);
            return Matrix4x4(T).init(
                f / aspect,
                0,
                0,
                0,
                0,
                f,
                0,
                0,
                0,
                0,
                far / (near - far),
                near * far / (near - far),
                0,
                0,
                -1,
                0,
            );
        }

        /// Generate the inverse perspective projection matrix.
        ///
        /// **Parameters**:
        /// - `A`: The perspective projection matrix.
        ///
        /// **Returns**:
        /// - The inverse perspective projection matrix.
        pub inline fn perspectiveInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1 / A.m[0].v[0],
                0,
                0,
                0,
                0,
                1 / A.m[1].v[1],
                0,
                0,
                0,
                0,
                0,
                A.m[2].v[3],
                0,
                0,
                1 / A.m[3].v[2],
                A.m[2].v[2] / A.m[3].v[2],
            );
        }

        /// Generate the orthographic projection matrix for a
        /// negative-one-to-one depth range. The projection matrix transforms
        /// from camera space to clip space.
        ///
        /// **Parameters**:
        /// - `left`: The left edge of the clipping plane.
        /// - `right`: The right edge of the clipping plane.
        /// - `bottom`: The bottom edge of the clipping plane.
        /// - `top`: The top edge of the clipping plane.
        /// - `near`: The distance to the near clipping plane.
        /// - `far`: The distance to the far clipping plane.
        ///
        /// **Returns**:
        /// - The orthographic projection matrix.
        pub inline fn orthographicNO(left: T, right: T, bottom: T, top: T, near: T, far: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                2 / (right - left),
                0,
                0,
                (right + left) / (left - right),
                0,
                2 / (top - bottom),
                0,
                (top + bottom) / (bottom - top),
                0,
                0,
                2 / (near - far),
                (near + far) / (near - far),
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the orthographic projection matrix for a  zero-to-one depth
        /// range. The projection matrix transforms from camera space to clip
        /// space.
        ///
        /// **Parameters**:
        /// - `left`: The left edge of the clipping plane.
        /// - `right`: The right edge of the clipping plane.
        /// - `bottom`: The bottom edge of the clipping plane.
        /// - `top`: The top edge of the clipping plane.
        /// - `near`: The distance to the near clipping plane.
        /// - `far`: The distance to the far clipping plane.
        ///
        /// **Returns**:
        /// - The orthographic projection matrix.
        pub inline fn orthographicZO(left: T, right: T, bottom: T, top: T, near: T, far: T) Matrix4x4(T) {
            return Matrix4x4(T).init(
                2 / (right - left),
                0,
                0,
                (right + left) / (left - right),
                0,
                2 / (top - bottom),
                0,
                (top + bottom) / (bottom - top),
                0,
                0,
                1 / (near - far),
                near / (near - far),
                0,
                0,
                0,
                1,
            );
        }

        /// Generate the inverse orthographic projection matrix.
        ///
        /// **Parameters**:
        /// - `A`: The orthographic projection matrix.
        ///
        /// **Returns**:
        /// - The inverse orthographic projection matrix.
        pub inline fn orthographicInverse(A: *const Matrix4x4(T)) Matrix4x4(T) {
            return Matrix4x4(T).init(
                1 / A.m[0].v[0],
                0,
                0,
                -A.m[3].v[0] / A.m[0].v[0],
                0,
                1 / A.m[1].v[1],
                0,
                -A.m[3].v[1] / A.m[1].v[1],
                0,
                0,
                1 / A.m[2].v[2],
                -A.m[3].v[2] / A.m[2].v[2],
                0,
                0,
                0,
                1,
            );
        }

        /// Creates the transformation matrix and applies it to the matrix
        /// in-place.
        ///
        /// **Parameters**:
        /// - `A`: The matrix to transform.
        /// - `f`: The transformation matrix generator function.
        /// - `args`: The function arguments.
        ///
        /// **Returns**:
        /// - The transformed matrix.
        pub inline fn applyTransform(A: *Matrix4x4(T), f: anytype, args: anytype) void {
            A.* = @call(.auto, f, args).mulMatrix4x4(A);
        }

        /// Compares two matrices and returns true if they are equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        ///
        /// **Returns**:
        /// - `true` if the matrices are equal, `false` otherwise.
        pub inline fn equal(A: *const Matrix4x4(T), B: *const Matrix4x4(T)) bool {
            return Vector4(T).equal(&A.m[0], &B.m[0]) and
                Vector4(T).equal(&A.m[1], &B.m[1]) and
                Vector4(T).equal(&A.m[2], &B.m[2]) and
                Vector4(T).equal(&A.m[3], &B.m[3]);
        }

        /// Compares two matrices and returns true if they are approximately
        /// equal.
        ///
        /// **Parameters**:
        /// - `A`: The first matrix.
        /// - `B`: The second matrix.
        /// - `tolerance`: The maximum difference between the components.
        ///
        /// **Returns**:
        /// - `true` if the matrices are approximately equal, `false` otherwise.
        pub inline fn approxEqual(A: *const Matrix4x4(T), B: *const Matrix4x4(T), tolerance: T) bool {
            return Vector4(T).approxEqual(&A.m[0], &B.m[0], tolerance) and
                Vector4(T).approxEqual(&A.m[1], &B.m[1], tolerance) and
                Vector4(T).approxEqual(&A.m[2], &B.m[2], tolerance) and
                Vector4(T).approxEqual(&A.m[3], &B.m[3], tolerance);
        }
    };
}

test "Matrix2x2.add" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x2(f32).init(
        5,
        6,
        7,
        8,
    );
    const R = Matrix2x2(f32).add(&A, &B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        6,
        8,
        10,
        12,
    ), 0.0001));
}

test "Matrix2x2.addScalar" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = Matrix2x2(f32).addScalar(&A, 5);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        6,
        7,
        8,
        9,
    ), 0.0001));
}

test "Matrix2x2.sub" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x2(f32).init(
        5,
        6,
        7,
        8,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        -4,
        -4,
        -4,
        -4,
    ), 0.0001));
}

test "Matrix2x2.subScalar" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.subScalar(5);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix2x2.negate" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        -1,
        -2,
        -3,
        -4,
    ), 0.0001));
}

test "Matrix2x2.mul" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x2(f32).init(
        5,
        6,
        7,
        8,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        5,
        12,
        21,
        32,
    ), 0.0001));
}

test "Matrix2x2.mulScalar" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.mulScalar(5);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        5,
        10,
        15,
        20,
    ), 0.0001));
}

test "Matrix2x2.mulVector2" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const v = Vector2(f32).init(5, 6);
    const R = A.mulVector2(&v);

    try std.testing.expect(R.approxEqual(&Vector2(f32).init(17, 39), 0.0001));
}

test "Matrix2x2.mulMatrix2x2" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x2(f32).init(
        5,
        6,
        7,
        8,
    );
    const R = A.mulMatrix2x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        19,
        22,
        43,
        50,
    ), 0.0001));
}

test "Matrix2x2.mulMatrix2x3" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x3(f32).init(
        5,
        6,
        7,
        8,
        9,
        10,
    );
    const R = A.mulMatrix2x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        21,
        24,
        27,
        47,
        54,
        61,
    ), 0.0001));
}

test "Matrix2x2.mulMatrix2x4" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x4(f32).init(
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mulMatrix2x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        23,
        26,
        29,
        32,
        51,
        58,
        65,
        72,
    ), 0.0001));
}

test "Matrix2x2.div" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const B = Matrix2x2(f32).init(
        5,
        6,
        7,
        8,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        0.2,
        0.3333,
        0.4285,
        0.5,
    ), 0.0001));
}

test "Matrix2x2.divScalar" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.divScalar(5);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        0.2,
        0.4,
        0.6,
        0.8,
    ), 0.0001));
}

test "Matrix2x2.transpose" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        1,
        3,
        2,
        4,
    ), 0.0001));
}

test "Matrix2x2.deteterminant" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.determinant();

    try std.testing.expectApproxEqRel(-2, R, 0.0001);
}

test "Matrix2x2.inverse" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = try A.inverse();

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        -2,
        1,
        1.5,
        -0.5,
    ), 0.0001));
}

test "Matrix2x2.trace" {
    const A = Matrix2x2(f32).init(
        1,
        2,
        3,
        4,
    );
    const R = A.trace();

    try std.testing.expectApproxEqRel(5, R, 0.0001);
}

test "Matrix2x3.add" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x3(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        8,
        10,
        12,
        14,
        16,
        18,
    ), 0.0001));
}

test "Matrix2x3.addScalar" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.addScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        8,
        9,
        10,
        11,
        12,
        13,
    ), 0.0001));
}

test "Matrix2x3.sub" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x3(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        -6,
        -6,
        -6,
        -6,
        -6,
        -6,
    ), 0.0001));
}

test "Matrix2x3.subScalar" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.subScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix2x3.negate" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
    ), 0.0001));
}

test "Matrix2x3.mul" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x3(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        7,
        16,
        27,
        40,
        55,
        72,
    ), 0.0001));
}

test "Matrix2x3.mulScalar" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.mulScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        7,
        14,
        21,
        28,
        35,
        42,
    ), 0.0001));
}

test "Matrix2x3.mulVector3" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const v = Vector3(f32).init(7, 8, 9);
    const R = A.mulVector3(&v);

    try std.testing.expect(R.approxEqual(&Vector2(f32).init(50, 122), 0.0001));
}

test "Matrix2x3.mulMatrix3x2" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x2(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mulMatrix3x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        58,
        64,
        139,
        154,
    ), 0.0001));
}

test "Matrix2x3.mulMatrix3x3" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x3(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
    );
    const R = A.mulMatrix3x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        66,
        72,
        78,
        156,
        171,
        186,
    ), 0.0001));
}

test "Matrix2x3.mulMatrix3x4" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x4(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.mulMatrix3x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        74,
        80,
        86,
        92,
        173,
        188,
        203,
        218,
    ), 0.0001));
}

test "Matrix2x3.div" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x3(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        0.1428,
        0.25,
        0.3333,
        0.4,
        0.4545,
        0.5,
    ), 0.0001));
}

test "Matrix2x3.divScalar" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.divScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        0.1428,
        0.2857,
        0.4285,
        0.5714,
        0.7142,
        0.8571,
    ), 0.0001));
}

test "Matrix2x3.transpose" {
    const A = Matrix2x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        1,
        4,
        2,
        5,
        3,
        6,
    ), 0.0001));
}

test "Matrix2x4.add" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x4(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        10,
        12,
        14,
        16,
        18,
        20,
        22,
        24,
    ), 0.0001));
}

test "Matrix2x4.addScalar" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.addScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
    ), 0.0001));
}

test "Matrix2x4.sub" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x4(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        -8,
        -8,
        -8,
        -8,
        -8,
        -8,
        -8,
        -8,
    ), 0.0001));
}

test "Matrix2x4.subScalar" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.subScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        -8,
        -7,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix2x4.negate" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
        -7,
        -8,
    ), 0.0001));
}

test "Matrix2x4.mul" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x4(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        9,
        20,
        33,
        48,
        65,
        84,
        105,
        128,
    ), 0.0001));
}

test "Matrix2x4.mulScalar" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.mulScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        9,
        18,
        27,
        36,
        45,
        54,
        63,
        72,
    ), 0.0001));
}

test "Matrix2x4.mulVector4" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const v = Vector4(f32).init(9, 10, 11, 12);
    const R = A.mulVector4(&v);

    try std.testing.expect(R.approxEqual(&Vector2(f32).init(110, 278), 0.0001));
}

test "Matrix2x4.mulMatrix4x2" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x2(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.mulMatrix4x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x2(f32).init(
        130,
        140,
        322,
        348,
    ), 0.0001));
}

test "Matrix2x4.mulMatrix4x3" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x3(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
    );
    const R = A.mulMatrix4x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        150,
        160,
        170,
        366,
        392,
        418,
    ), 0.0001));
}

test "Matrix2x4.mulMatrix4x4" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x4(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.mulMatrix4x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        170,
        180,
        190,
        200,
        410,
        436,
        462,
        488,
    ), 0.0001));
}

test "Matrix2x4.div" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x4(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        0.1111,
        0.2,
        0.2727,
        0.3333,
        0.3846,
        0.4285,
        0.4666,
        0.5,
    ), 0.0001));
}

test "Matrix2x4.divScalar" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.divScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        0.1111,
        0.2222,
        0.3333,
        0.4444,
        0.5555,
        0.6666,
        0.7777,
        0.8888,
    ), 0.0001));
}

test "Matrix2x4.transpose" {
    const A = Matrix2x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        1,
        5,
        2,
        6,
        3,
        7,
        4,
        8,
    ), 0.0001));
}

test "Matrix3x2.add" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x2(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        8,
        10,
        12,
        14,
        16,
        18,
    ), 0.0001));
}

test "Matrix3x2.addScalar" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.addScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        8,
        9,
        10,
        11,
        12,
        13,
    ), 0.0001));
}

test "Matrix3x2.sub" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x2(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        -6,
        -6,
        -6,
        -6,
        -6,
        -6,
    ), 0.0001));
}

test "Matrix3x2.subScalar" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.subScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix3x2.negate" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
    ), 0.0001));
}

test "Matrix3x2.mul" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x2(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        7,
        16,
        27,
        40,
        55,
        72,
    ), 0.0001));
}

test "Matrix3x2.mulScalar" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.mulScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        7,
        14,
        21,
        28,
        35,
        42,
    ), 0.0001));
}

test "Matrix3x2.mulVector2" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const v = Vector2(f32).init(7, 8);
    const R = A.mulVector2(&v);

    try std.testing.expect(R.approxEqual(&Vector3(f32).init(23, 53, 83), 0.0001));
}

test "Matrix3x2.mulMatrix2x2" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x2(f32).init(
        7,
        8,
        9,
        10,
    );
    const R = A.mulMatrix2x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        25,
        28,
        57,
        64,
        89,
        100,
    ), 0.0001));
}

test "Matrix3x2.mulMatrix2x3" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x3(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mulMatrix2x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        27,
        30,
        33,
        61,
        68,
        75,
        95,
        106,
        117,
    ), 0.0001));
}

test "Matrix3x2.mulMatrix2x4" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix2x4(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
    );
    const R = A.mulMatrix2x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        29,
        32,
        35,
        38,
        65,
        72,
        79,
        86,
        101,
        112,
        123,
        134,
    ), 0.0001));
}

test "Matrix3x2.div" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const B = Matrix3x2(f32).init(
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        0.1428,
        0.25,
        0.3333,
        0.4,
        0.4545,
        0.5,
    ), 0.0001));
}

test "Matrix3x2.divScalar" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.divScalar(7);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        0.1428,
        0.2857,
        0.4285,
        0.5714,
        0.7142,
        0.8571,
    ), 0.0001));
}

test "Matrix3x2.transpose" {
    const A = Matrix3x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix2x3(f32).init(
        1,
        3,
        5,
        2,
        4,
        6,
    ), 0.0001));
}

test "Matrix3x3.add" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x3(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        11,
        13,
        15,
        17,
        19,
        21,
        23,
        25,
        27,
    ), 0.0001));
}

test "Matrix3x3.addScalar" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.addScalar(10);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
    ), 0.0001));
}

test "Matrix3x3.sub" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x3(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        -9,
        -9,
        -9,
        -9,
        -9,
        -9,
        -9,
        -9,
        -9,
    ), 0.0001));
}

test "Matrix3x3.subScalar" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.subScalar(10);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        -9,
        -8,
        -7,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix3x3.negate" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
        -7,
        -8,
        -9,
    ), 0.0001));
}

test "Matrix3x3.mul" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x3(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        10,
        22,
        36,
        52,
        70,
        90,
        112,
        136,
        162,
    ), 0.0001));
}

test "Matrix3x3.mulScalar" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.mulScalar(10);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        80,
        90,
    ), 0.0001));
}

test "Matrix3x3.mulVector2" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const v = Vector2(f32).init(10, 11);
    const R = A.mulVector2(&v);

    try std.testing.expect(R.approxEqual(&Vector2(f32).init(35, 101), 0.0001));
}

test "Matrix3x3.mulVector3" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const v = Vector3(f32).init(10, 11, 12);
    const R = A.mulVector3(&v);

    try std.testing.expect(R.approxEqual(&Vector3(f32).init(68, 167, 266), 0.0001));
}

test "Matrix3x3.mulMatrix3x2" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x2(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
    );
    const R = A.mulMatrix3x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        76,
        82,
        184,
        199,
        292,
        316,
    ), 0.0001));
}

test "Matrix3x3.mulMatrix3x3" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x3(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.mulMatrix3x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        84,
        90,
        96,
        201,
        216,
        231,
        318,
        342,
        366,
    ), 0.0001));
}

test "Matrix3x3.mulMatrix3x4" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x4(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
    );
    const R = A.mulMatrix3x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        92,
        98,
        104,
        110,
        218,
        233,
        248,
        263,
        344,
        368,
        392,
        416,
    ), 0.0001));
}

test "Matrix3x3.div" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const B = Matrix3x3(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        0.1,
        0.1818,
        0.25,
        0.3076,
        0.3571,
        0.4,
        0.4375,
        0.4705,
        0.5,
    ), 0.0001));
}

test "Matrix3x3.divScalar" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.divScalar(10);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
    ), 0.0001));
}

test "Matrix3x3.transpose" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        1,
        4,
        7,
        2,
        5,
        8,
        3,
        6,
        9,
    ), 0.0001));
}

test "Matrix3x3.deteterminant" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.determinant();

    try std.testing.expectApproxEqRel(0, R, 0.0001);
}

test "Matrix3x3.inverse" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        0,
        1,
        4,
        5,
        6,
        0,
    );
    const R = try A.inverse();

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        -24,
        18,
        5,
        20,
        -15,
        -4,
        -5,
        4,
        1,
    ), 0.0001));
}

test "Matrix3x3.trace" {
    const A = Matrix3x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
    );
    const R = A.trace();

    try std.testing.expectApproxEqRel(15, R, 0.0001);
}

test "Matrix3x3.scale" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).scale(6, 7);
    const Ainv = A.scaleInverse();
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(6, 14), 0.0001));

    R1 = Ainv.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(18, 28, 1), 0.0001));

    R2 = Ainv.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.shearX" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).shearX(5);
    const Ainv = A.shearInverse();
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(11, 2), 0.0001));

    R1 = Ainv.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(23, 4, 1), 0.0001));

    R2 = Ainv.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.shearY" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).shearY(5);
    const Ainv = A.shearInverse();
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(1, 7), 0.0001));

    R1 = Ainv.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(3, 19, 1), 0.0001));

    R2 = Ainv.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.rotate" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).rotate(3.14159 / 2.0);
    const Ainv = A.rotateInverse();
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(2, -1), 0.0001));

    R1 = Ainv.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(4, -3, 1), 0.0001));

    R2 = Ainv.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.reflect" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).reflect(&Vector2(f32).init(1, 0).normalize());
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(-1, 2), 0.0001));

    R1 = A.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(-3, 4, 1), 0.0001));

    R2 = A.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.reflectPoint" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).reflectPoint(&Vector2(f32).init(0, 1).normalize(), &Vector2(f32).init(1, 1));
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(1, 0), 0.0001));

    R1 = A.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(3, -2, 1), 0.0001));

    R2 = A.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.translate" {
    const v = Vector2(f32).init(1, 2);
    const w = Vector3(f32).init(3, 4, 1);
    const A = Matrix3x3(f32).translate(5, 6);
    const Ainv = A.translateInverse();
    var R1 = A.mulVector2(&v);

    try std.testing.expect(R1.approxEqual(&Vector2(f32).init(6, 8), 0.0001));

    R1 = Ainv.mulVector2(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector3(&w);

    try std.testing.expect(R2.approxEqual(&Vector3(f32).init(8, 10, 1), 0.0001));

    R2 = Ainv.mulVector3(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix3x3.applyTransform" {
    var A = Matrix3x3(f32).identity;
    A.applyTransform(Matrix3x3(f32).scale, .{ 2.0, 3.0 });

    try std.testing.expect(A.approxEqual(&Matrix3x3(f32).scale(2.0, 3.0), 0.0001));
}

test "Matrix3x4.add" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x4(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        14,
        16,
        18,
        20,
        22,
        24,
        26,
        28,
        30,
        32,
        34,
        36,
    ), 0.0001));
}

test "Matrix3x4.addScalar" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.addScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ), 0.0001));
}

test "Matrix3x4.sub" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x4(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
    ), 0.0001));
}

test "Matrix3x4.subScalar" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.subScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        -12,
        -11,
        -10,
        -9,
        -8,
        -7,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix3x4.negate" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
        -7,
        -8,
        -9,
        -10,
        -11,
        -12,
    ), 0.0001));
}

test "Matrix3x4.mul" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x4(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        13,
        28,
        45,
        64,
        85,
        108,
        133,
        160,
        189,
        220,
        253,
        288,
    ), 0.0001));
}

test "Matrix3x4.mulScalar" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mulScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        13,
        26,
        39,
        52,
        65,
        78,
        91,
        104,
        117,
        130,
        143,
        156,
    ), 0.0001));
}

test "Matrix3x4.mulVector4" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const v = Vector4(f32).init(13, 14, 15, 16);
    const R = A.mulVector4(&v);

    try std.testing.expect(R.approxEqual(&Vector3(f32).init(150, 382, 614), 0.0001));
}

test "Matrix3x4.mulMatrix4x2" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x2(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
    );
    const R = A.mulMatrix4x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x2(f32).init(
        170,
        180,
        426,
        452,
        682,
        724,
    ), 0.0001));
}

test "Matrix3x4.mulMatrix4x3" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x3(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.mulMatrix4x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x3(f32).init(
        190,
        200,
        210,
        470,
        496,
        522,
        750,
        792,
        834,
    ), 0.0001));
}

test "Matrix3x4.mulMatrix4x4" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x4(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
    );
    const R = A.mulMatrix4x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        210,
        220,
        230,
        240,
        514,
        540,
        566,
        592,
        818,
        860,
        902,
        944,
    ), 0.0001));
}

test "Matrix3x4.div" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x4(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        0.0769,
        0.1428,
        0.2,
        0.25,
        0.2941,
        0.3333,
        0.3684,
        0.4,
        0.4285,
        0.4545,
        0.4782,
        0.5,
    ), 0.0001));
}

test "Matrix3x4.divScalar" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.divScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        0.0769,
        0.1538,
        0.2307,
        0.3076,
        0.3846,
        0.4615,
        0.5384,
        0.6153,
        0.6923,
        0.7692,
        0.8461,
        0.923,
    ), 0.0001));
}

test "Matrix3x4.transpose" {
    const A = Matrix3x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        1,
        5,
        9,
        2,
        6,
        10,
        3,
        7,
        11,
        4,
        8,
        12,
    ), 0.0001));
}

test "Matrix4x2.add" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x2(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        10,
        12,
        14,
        16,
        18,
        20,
        22,
        24,
    ), 0.0001));
}

test "Matrix4x2.addScalar" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.addScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
    ), 0.0001));
}

test "Matrix4x2.sub" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x2(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        -8,
        -8,
        -8,
        -8,
        -8,
        -8,
        -8,
        -8,
    ), 0.0001));
}

test "Matrix4x2.subScalar" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.subScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        -8,
        -7,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix4x2.negate" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
        -7,
        -8,
    ), 0.0001));
}

test "Matrix4x2.mul" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x2(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        9,
        20,
        33,
        48,
        65,
        84,
        105,
        128,
    ), 0.0001));
}

test "Matrix4x2.mulScalar" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.mulScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        9,
        18,
        27,
        36,
        45,
        54,
        63,
        72,
    ), 0.0001));
}

test "Matrix4x2.mulVector2" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const v = Vector2(f32).init(9, 10);
    const R = A.mulVector2(&v);

    try std.testing.expect(R.approxEqual(&Vector4(f32).init(29, 67, 105, 143), 0.0001));
}

test "Matrix4x2.mulMatrix2x2" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x2(f32).init(
        9,
        10,
        11,
        12,
    );
    const R = A.mulMatrix2x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        31,
        34,
        71,
        78,
        111,
        122,
        151,
        166,
    ), 0.0001));
}

test "Matrix4x2.mulMatrix2x3" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x3(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
    );
    const R = A.mulMatrix2x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        33,
        36,
        39,
        75,
        82,
        89,
        117,
        128,
        139,
        159,
        174,
        189,
    ), 0.0001));
}

test "Matrix4x2.mulMatrix2x4" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix2x4(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.mulMatrix2x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        35,
        38,
        41,
        44,
        79,
        86,
        93,
        100,
        123,
        134,
        145,
        156,
        167,
        182,
        197,
        212,
    ), 0.0001));
}

test "Matrix4x2.div" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const B = Matrix4x2(f32).init(
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        0.1111,
        0.2,
        0.2727,
        0.3333,
        0.3846,
        0.4285,
        0.4666,
        0.5,
    ), 0.0001));
}

test "Matrix4x2.divScalar" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.divScalar(9);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        0.1111,
        0.2222,
        0.3333,
        0.4444,
        0.5555,
        0.6666,
        0.7777,
        0.8888,
    ), 0.0001));
}

test "Matrix4x2.transpose" {
    const A = Matrix4x2(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix2x4(f32).init(
        1,
        3,
        5,
        7,
        2,
        4,
        6,
        8,
    ), 0.0001));
}

test "Matrix4x3.add" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x3(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        14,
        16,
        18,
        20,
        22,
        24,
        26,
        28,
        30,
        32,
        34,
        36,
    ), 0.0001));
}

test "Matrix4x3.addScalar" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.addScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ), 0.0001));
}

test "Matrix4x3.sub" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x3(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
        -12,
    ), 0.0001));
}

test "Matrix4x3.subScalar" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.subScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        -12,
        -11,
        -10,
        -9,
        -8,
        -7,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix4x3.negate" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
        -7,
        -8,
        -9,
        -10,
        -11,
        -12,
    ), 0.0001));
}

test "Matrix4x3.mul" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x3(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        13,
        28,
        45,
        64,
        85,
        108,
        133,
        160,
        189,
        220,
        253,
        288,
    ), 0.0001));
}

test "Matrix4x3.mulScalar" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.mulScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        13,
        26,
        39,
        52,
        65,
        78,
        91,
        104,
        117,
        130,
        143,
        156,
    ), 0.0001));
}

test "Matrix4x3.mulVector3" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const v = Vector3(f32).init(13, 14, 15);
    const R = A.mulVector3(&v);

    try std.testing.expect(R.approxEqual(&Vector4(f32).init(86, 212, 338, 464), 0.0001));
}

test "Matrix4x3.mulMatrix3x2" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x2(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
    );
    const R = A.mulMatrix3x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        94,
        100,
        229,
        244,
        364,
        388,
        499,
        532,
    ), 0.0001));
}

test "Matrix4x3.mulMatrix3x3" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x3(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
    );
    const R = A.mulMatrix3x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        102,
        108,
        114,
        246,
        261,
        276,
        390,
        414,
        438,
        534,
        567,
        600,
    ), 0.0001));
}

test "Matrix4x3.mulMatrix3x4" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix3x4(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.mulMatrix3x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        110,
        116,
        122,
        128,
        263,
        278,
        293,
        308,
        416,
        440,
        464,
        488,
        569,
        602,
        635,
        668,
    ), 0.0001));
}

test "Matrix4x3.div" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const B = Matrix4x3(f32).init(
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        0.0769,
        0.1428,
        0.2,
        0.25,
        0.2941,
        0.3333,
        0.3684,
        0.4,
        0.4285,
        0.4545,
        0.4782,
        0.5,
    ), 0.0001));
}

test "Matrix4x3.divScalar" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.divScalar(13);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        0.0769,
        0.1538,
        0.2307,
        0.3076,
        0.3846,
        0.4615,
        0.5384,
        0.6153,
        0.6923,
        0.7692,
        0.8461,
        0.923,
    ), 0.0001));
}

test "Matrix4x3.transpose" {
    const A = Matrix4x3(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix3x4(f32).init(
        1,
        4,
        7,
        10,
        2,
        5,
        8,
        11,
        3,
        6,
        9,
        12,
    ), 0.0001));
}

test "Matrix4x4.add" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x4(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
    );
    const R = A.add(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        18,
        20,
        22,
        24,
        26,
        28,
        30,
        32,
        34,
        36,
        38,
        40,
        42,
        44,
        46,
        48,
    ), 0.0001));
}

test "Matrix4x4.addScalar" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.addScalar(17);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
    ), 0.0001));
}

test "Matrix4x4.sub" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x4(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
    );
    const R = A.sub(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
        -16,
    ), 0.0001));
}

test "Matrix4x4.subScalar" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.subScalar(17);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        -16,
        -15,
        -14,
        -13,
        -12,
        -11,
        -10,
        -9,
        -8,
        -7,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
    ), 0.0001));
}

test "Matrix4x4.negate" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.negate();

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        -1,
        -2,
        -3,
        -4,
        -5,
        -6,
        -7,
        -8,
        -9,
        -10,
        -11,
        -12,
        -13,
        -14,
        -15,
        -16,
    ), 0.0001));
}

test "Matrix4x4.mul" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x4(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
    );
    const R = A.mul(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        17,
        36,
        57,
        80,
        105,
        132,
        161,
        192,
        225,
        260,
        297,
        336,
        377,
        420,
        465,
        512,
    ), 0.0001));
}

test "Matrix4x4.mulScalar" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.mulScalar(17);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        17,
        34,
        51,
        68,
        85,
        102,
        119,
        136,
        153,
        170,
        187,
        204,
        221,
        238,
        255,
        272,
    ), 0.0001));
}

test "Matrix4x4.mulVector3" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const v = Vector3(f32).init(17, 18, 19);
    const R = A.mulVector3(&v);

    try std.testing.expect(R.approxEqual(&Vector3(f32).init(114, 334, 554), 0.0001));
}

test "Matrix4x4.mulVector4" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const v = Vector4(f32).init(17, 18, 19, 20);
    const R = A.mulVector4(&v);

    try std.testing.expect(R.approxEqual(&Vector4(f32).init(190, 486, 782, 1078), 0.0001));
}

test "Matrix4x4.mulMatrix4x2" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x2(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    );
    const R = A.mulMatrix4x2(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x2(f32).init(
        210,
        220,
        530,
        556,
        850,
        892,
        1170,
        1228,
    ), 0.0001));
}

test "Matrix4x4.mulMatrix4x3" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x3(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
    );
    const R = A.mulMatrix4x3(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x3(f32).init(
        230,
        240,
        250,
        574,
        600,
        626,
        918,
        960,
        1002,
        1262,
        1320,
        1378,
    ), 0.0001));
}

test "Matrix4x4.mulMatrix4x4" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x4(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
    );
    const R = A.mulMatrix4x4(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        250,
        260,
        270,
        280,
        618,
        644,
        670,
        696,
        986,
        1028,
        1070,
        1112,
        1354,
        1412,
        1470,
        1528,
    ), 0.0001));
}

test "Matrix4x4.div" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const B = Matrix4x4(f32).init(
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
    );
    const R = A.div(&B);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        0.0588,
        0.1111,
        0.1578,
        0.2,
        0.2381,
        0.2727,
        0.3043,
        0.3333,
        0.36,
        0.3846,
        0.4074,
        0.4285,
        0.4482,
        0.4666,
        0.4838,
        0.5,
    ), 0.0001));
}

test "Matrix4x4.divScalar" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.divScalar(17);

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        0.0588,
        0.1176,
        0.1764,
        0.2352,
        0.2941,
        0.3529,
        0.4117,
        0.4705,
        0.5294,
        0.5882,
        0.647,
        0.7058,
        0.7647,
        0.8235,
        0.8823,
        0.9411,
    ), 0.0001));
}

test "Matrix4x4.transpose" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.transpose();

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        1,
        5,
        9,
        13,
        2,
        6,
        10,
        14,
        3,
        7,
        11,
        15,
        4,
        8,
        12,
        16,
    ), 0.0001));
}

test "Matrix4x4.deteterminant" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.determinant();

    try std.testing.expectApproxEqRel(0, R, 0.0001);
}

test "Matrix4x4.inverse" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        0,
        1,
        5,
        6,
        7,
        8,
        0,
        9,
        10,
        11,
        12,
        0,
    );
    const R = try A.inverse();

    try std.testing.expect(R.approxEqual(&Matrix4x4(f32).init(
        -2.8918,
        1.3963,
        0.3543,
        0.1411,
        2.8648,
        -1.4954,
        -0.2762,
        -0.0930,
        -0.2162,
        0.2072,
        -0.0420,
        0.0510,
        -0.2972,
        0.2432,
        0.0810,
        -0.0270,
    ), 0.0001));
}

test "Matrix4x4.trace" {
    const A = Matrix4x4(f32).init(
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
    );
    const R = A.trace();

    try std.testing.expectApproxEqRel(34, R, 0.0001);
}

test "Matrix4x4.scale" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).scale(7, 8, 9);
    const Ainv = A.scaleInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(7, 16, 27), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(28, 40, 54, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.shearX" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).shearX(7, 8);
    const Ainv = A.shearInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(39, 2, 3), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(87, 5, 6, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.shearY" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).shearY(7, 8);
    const Ainv = A.shearInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(1, 33, 3), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(4, 81, 6, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.shearZ" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).shearZ(7, 8);
    const Ainv = A.shearInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(1, 2, 26), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(4, 5, 74, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.rotate" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).rotate(3.14159 / 2.0, &Vector3(f32).init(0, 1, 0).normalize());
    const Ainv = A.rotateInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(3, 2, -1), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(6, 5, -4, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.rotateX" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).rotateX(3.14159 / 2.0);
    const Ainv = A.rotateInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(1, -3, 2), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(4, -6, 5, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.rotateY" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).rotateY(3.14159 / 2.0);
    const Ainv = A.rotateInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(3, 2, -1), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(6, 5, -4, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.rotateZ" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).rotateZ(3.14159 / 2.0);
    const Ainv = A.rotateInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(-2, 1, 3), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(-5, 4, 6, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.reflect" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const N = Vector3(f32).init(0, 1, 0);
    const A = Matrix4x4(f32).reflect(&N.normalize());
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(1, -2, 3), 0.0001));

    R1 = A.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(4, -5, 6, 1), 0.0001));

    R2 = A.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.reflectPoint" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const N = Vector3(f32).init(0, 1, 0);
    const p = Vector3(f32).init(0, 1, 0);
    const A = Matrix4x4(f32).reflectPoint(&N.normalize(), &p);
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(1, 0, 3), 0.0001));

    R1 = A.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(4, -3, 6, 1), 0.0001));

    R2 = A.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.translate" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).translate(7, 8, 9);
    const Ainv = A.translateInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(8, 10, 12), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(11, 13, 15, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.lookAt" {
    const v = Vector3(f32).init(1, 2, 3);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const eye = Vector3(f32).init(0, 0, 0);
    const center = Vector3(f32).init(0, 1, 1);
    const up = Vector3(f32).init(0, 1, 0);
    const A = Matrix4x4(f32).lookAt(&eye, &center, &up);
    const Ainv = A.lookAtInverse();
    var R1 = A.mulVector3(&v);

    try std.testing.expect(R1.approxEqual(&Vector3(f32).init(-1, -0.7071, -3.5355), 0.0001));

    R1 = Ainv.mulVector3(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(-4, -0.7071, -7.7781, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.perspectiveNO" {
    const v = Vector4(f32).init(1, 2, 3, 1);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).perspectiveNO(std.math.degreesToRadians(90.0), 16.0 / 9.0, 0.1, 100.0);
    const Ainv = A.perspectiveInverse();
    var R1 = A.mulVector4(&v);

    try std.testing.expect(R1.approxEqual(&Vector4(f32).init(0.5625, 2, -3.2062, -3), 0.0001));

    R1 = Ainv.mulVector4(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(2.25, 5, -6.2122, -6), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.perspectiveZO" {
    const v = Vector4(f32).init(1, 2, 3, 1);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).perspectiveZO(std.math.degreesToRadians(90.0), 16.0 / 9.0, 0.1, 100.0);
    const Ainv = A.perspectiveInverse();
    var R1 = A.mulVector4(&v);

    try std.testing.expect(R1.approxEqual(&Vector4(f32).init(0.5625, 2, -3.1031, -3), 0.0001));

    R1 = Ainv.mulVector4(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(2.25, 5, -6.1061, -6), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.orthographicNO" {
    const v = Vector4(f32).init(1, 2, 3, 1);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).orthographicNO(-10.0, 10.0, -10.0, 10.0, 0.1, 100.0);
    const Ainv = A.orthographicInverse();
    var R1 = A.mulVector4(&v);

    try std.testing.expect(R1.approxEqual(&Vector4(f32).init(0.1, 0.2, -1.0620, 1), 0.0001));

    R1 = Ainv.mulVector4(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(0.4, 0.5, -1.1221, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.orthographicZO" {
    const v = Vector4(f32).init(1, 2, 3, 1);
    const w = Vector4(f32).init(4, 5, 6, 1);
    const A = Matrix4x4(f32).orthographicZO(-10.0, 10.0, -10.0, 10.0, 0.1, 100.0);
    const Ainv = A.orthographicInverse();
    var R1 = A.mulVector4(&v);

    try std.testing.expect(R1.approxEqual(&Vector4(f32).init(0.1, 0.2, -0.0310, 1), 0.0001));

    R1 = Ainv.mulVector4(&R1);

    try std.testing.expect(R1.approxEqual(&v, 0.0001));

    var R2 = A.mulVector4(&w);

    try std.testing.expect(R2.approxEqual(&Vector4(f32).init(0.4, 0.5, -0.0610, 1), 0.0001));

    R2 = Ainv.mulVector4(&R2);

    try std.testing.expect(R2.approxEqual(&w, 0.0001));
}

test "Matrix4x4.applyTransform" {
    var A = Matrix4x4(f32).identity;
    A.applyTransform(Matrix4x4(f32).scale, .{ 2.0, 3.0, 4.0 });

    try std.testing.expect(A.approxEqual(&Matrix4x4(f32).scale(2.0, 3.0, 4.0), 0.0001));
}
