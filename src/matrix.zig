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
        /// **Returns**:
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
        pub inline fn initColumn(c0: *const Vector2(T), c1: *const Vector2(T), c2: *const Vector2(T)) Matrix2x3(T) {
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
        pub inline fn mulVector3(A: *const Matrix2x3(T), v: Vector3(T)) Vector2(T) {
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
        pub inline fn initColumn(c0: *const Vector2(T), c1: *const Vector2(T), c2: *const Vector2(T), c3: *const Vector2(T)) Matrix2x4(T) {
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
        pub inline fn mulVector4(A: *const Matrix2x4(T), v: Vector4(T)) Vector2(T) {
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
        pub inline fn initColumn(c0: *const Vector3(T), c1: *const Vector3(T)) Matrix3x2(T) {
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
        pub inline fn mulVector2(A: *const Matrix3x2(T), v: Vector2(T)) Vector3(T) {
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
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[0].v[1] * B.m[0].v[2],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[0].v[2] * B.m[0].v[2],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[0].v[3] * B.m[0].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[0].v[1] * B.m[1].v[2],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[0].v[2] * B.m[1].v[2],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[0].v[3] * B.m[1].v[2],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[0].v[1] * B.m[2].v[2],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[0].v[2] * B.m[2].v[2],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[0].v[3] * B.m[2].v[2],
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
                    A.m[0].v[0] * B.m[0].v[0] + A.m[1].v[0] * B.m[0].v[1] + A.m[0].v[1] * B.m[0].v[2] + A.m[1].v[1] * B.m[0].v[3],
                    A.m[0].v[1] * B.m[0].v[0] + A.m[1].v[1] * B.m[0].v[1] + A.m[0].v[2] * B.m[0].v[2] + A.m[1].v[2] * B.m[0].v[3],
                    A.m[0].v[2] * B.m[0].v[0] + A.m[1].v[2] * B.m[0].v[1] + A.m[0].v[3] * B.m[0].v[2] + A.m[1].v[3] * B.m[0].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[1].v[0] + A.m[1].v[0] * B.m[1].v[1] + A.m[0].v[1] * B.m[1].v[2] + A.m[1].v[1] * B.m[1].v[3],
                    A.m[0].v[1] * B.m[1].v[0] + A.m[1].v[1] * B.m[1].v[1] + A.m[0].v[2] * B.m[1].v[2] + A.m[1].v[2] * B.m[1].v[3],
                    A.m[0].v[2] * B.m[1].v[0] + A.m[1].v[2] * B.m[1].v[1] + A.m[0].v[3] * B.m[1].v[2] + A.m[1].v[3] * B.m[1].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[2].v[0] + A.m[1].v[0] * B.m[2].v[1] + A.m[0].v[1] * B.m[2].v[2] + A.m[1].v[1] * B.m[2].v[3],
                    A.m[0].v[1] * B.m[2].v[0] + A.m[1].v[1] * B.m[2].v[1] + A.m[0].v[2] * B.m[2].v[2] + A.m[1].v[2] * B.m[2].v[3],
                    A.m[0].v[2] * B.m[2].v[0] + A.m[1].v[2] * B.m[2].v[1] + A.m[0].v[3] * B.m[2].v[2] + A.m[1].v[3] * B.m[2].v[3],
                ),
                Vector3(T).init(
                    A.m[0].v[0] * B.m[3].v[0] + A.m[1].v[0] * B.m[3].v[1] + A.m[0].v[1] * B.m[3].v[2] + A.m[1].v[1] * B.m[3].v[3],
                    A.m[0].v[1] * B.m[3].v[0] + A.m[1].v[1] * B.m[3].v[1] + A.m[0].v[2] * B.m[3].v[2] + A.m[1].v[2] * B.m[3].v[3],
                    A.m[0].v[2] * B.m[3].v[0] + A.m[1].v[2] * B.m[3].v[1] + A.m[0].v[3] * B.m[3].v[2] + A.m[1].v[3] * B.m[3].v[3],
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
        pub inline fn initColumn(c0: *const Vector3(T), c1: *const Vector3(T), c2: *const Vector3(T)) Matrix3x3(T) {
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
        pub const eye = Matrix3x3(T).init(1, 0, 0, 0, 1, 0, 0, 0, 1);

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
        pub inline fn mulVector2(A: *const Matrix3x3(T), v: Vector2(T)) Vector2(T) {
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
        pub inline fn mulVector3(A: *const Matrix3x3(T), v: Vector3(T)) Vector3(T) {
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
        pub inline fn inverse(A: *const Matrix3x3(T)) Matrix3x3(T) {
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

        /// Generate the shearing matrix.
        ///
        /// **Parameters**:
        /// - `x`: The shearing factor along the x-axis.
        /// - `y`: The shearing factor along the y-axis.
        ///
        /// **Returns**:
        /// - The shearing matrix.
        pub inline fn shear(x: T, y: T) Matrix3x3(T) {
            return Matrix3x3(T).init(1, x, 0, y, 1, 0, 0, 0, 1);
        }

        /// Generate the inverse shearing matrix.
        ///
        /// **Parameters**:
        /// - `A`: The shearing matrix.
        ///
        /// **Returns**:
        /// - The inverse shearing matrix.
        pub inline fn shearInverse(A: *const Matrix3x3(T)) Matrix3x3(T) {
            return Matrix3x3(T).init(1, -A.m[0].v[1], 0, -A.m[1].v[0], 1, 0, 0, 0, 1);
        }

        /// Generate the left-handed rotation matrix.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotateLH(angle: T) Matrix3x3(T) {
            const c = @cos(angle);
            const s = @sin(angle);
            return Matrix3x3(T).init(c, -s, 0, s, c, 0, 0, 0, 1);
        }

        /// Generate the right-handed rotation matrix.
        ///
        /// **Parameters**:
        /// - `angle`: The rotation angle in radians.
        ///
        /// **Returns**:
        /// - The rotation matrix.
        pub inline fn rotateRH(angle: T) Matrix3x3(T) {
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
            return Matrix3x3(T).init(A.m[0].v[0], A.m[1].v[0], 0, A.m[0].v[1], A.m[1].v[1], 0, 0, 0, 1);
        }

        /// Generate the reflection matrix across a normal vector.
        ///
        /// **Parameters**:
        /// - `n`: The normal vector. Must be normalized.
        ///
        /// **Returns**:
        /// - The reflection matrix.
        pub inline fn reflect(normal: Vector2(T)) Matrix3x3(T) {
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
