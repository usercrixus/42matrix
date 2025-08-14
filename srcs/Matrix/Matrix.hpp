#pragma once

#include "VectorAlgebra.hpp"
#include <complex>

template <typename T>
class Matrix : private std::vector<VectorAlgebra<T>>
{
private:
    /**
     * @brief Check if the matrix is valid.
     * @return true if all rows have the same number of columns, false otherwise.
     */
    bool isValidate() const;

    /**
     * @brief Check if the matrix is square.
     * @return true if number of rows == number of columns, false otherwise.
     */
    bool isSquare() const;

public:
    using std::vector<VectorAlgebra<T>>::operator[];
    using std::vector<VectorAlgebra<T>>::begin;
    using std::vector<VectorAlgebra<T>>::end;

    /**
     * @brief Create a Matrix from a list of rows.
     * @param rows Vector of VectorAlgebra<T> representing the matrix rows.
     * @throws std::invalid_argument if rows are not all the same size.
     */
    static Matrix<T> from(const std::vector<VectorAlgebra<T>> &rows);

    /**
     * @brief Interpolate between two matrices.
     * @param m1 First matrix.
     * @param m2 Second matrix.
     * @param ratio Interpolation factor in [0, 1].
     * @return The interpolated matrix.
     * @throws std::invalid_argument if matrix dimensions differ.
     */
    static Matrix<T> linearInterpolation(const Matrix<T> &m1, const Matrix<T> &m2, float ratio);

    /**
     * @brief Create a 4x4 perspective projection matrix.
     * @param fov Field of view in degrees.
     * @param ratio Aspect ratio (width / height).
     * @param near Near clipping plane distance.
     * @param far Far clipping plane distance.
     * @return Projection matrix.
     */
    static Matrix<T> projection(float fov, float ratio, float near, float far);

    /**
     * @brief Get the transposed version of the matrix.
     * @return Transposed matrix.
     */
    Matrix<T> transpose() const;

    /**
     * @brief Convert matrix to row echelon form (Gauss-Jordan elimination).
     * @return Row echelon form of the matrix.
     */
    Matrix<T> rowEchelon() const;

    /**
     * @brief Row echelon form with pivot normalization and swap count tracking.
     * @param swapCount Output parameter for number of row swaps performed.
     * @return Row echelon normalized matrix.
     */
    Matrix<T> rowEchelonNormalize(int &swapCount) const;

    /**
     * @brief Compute the inverse of a square matrix.
     * @return Inverted matrix.
     * @throws std::logic_error if the matrix is not square or non-invertible.
     */
    Matrix<T> invert() const;

    /**
     * @brief Generate an identity matrix of the same size as this matrix.
     * @return Identity matrix.
     */
    Matrix<T> getIdentityMatrix() const;

    /**
     * @brief Compute the rank of the matrix.
     * @return Rank as an unsigned integer.
     */
    unsigned int rank() const;

    /**
     * @brief Compute the determinant of the matrix.
     * @return Determinant value.
     * @throws std::invalid_argument if the matrix is not square.
     */
    T determinant() const;

    /**
     * @brief Compute the trace of a square matrix (sum of diagonal elements).
     * @return Trace value.
     * @throws std::logic_error if the matrix is not square or empty.
     */
    T trace();

    /**
     * @brief Multiply this matrix by another matrix.
     * @param other Matrix to multiply with.
     * @return Product matrix.
     * @throws std::invalid_argument if dimensions are incompatible.
     */
    Matrix<T> operator*(const Matrix<T> &other) const;

    /**
     * @brief Multiply this matrix by a scalar.
     * @tparam Scalar Numeric type.
     * @param scalar Value to multiply with.
     * @return Scaled matrix.
     */
    template <typename Scalar>
    Matrix<T> operator*(const Scalar scalar) const;

    /**
     * @brief Add two matrices element-wise.
     * @param other Matrix to add.
     * @return Sum matrix.
     * @throws std::invalid_argument if dimensions differ.
     */
    Matrix<T> operator+(const Matrix<T> &other) const;

    /**
     * @brief Subtract another matrix from this matrix element-wise.
     * @param other Matrix to subtract.
     * @return Difference matrix.
     * @throws std::invalid_argument if dimensions differ.
     */
    Matrix<T> operator-(const Matrix<T> &other) const;

    /**
     * @brief Multiply this matrix by a vector.
     * @param vec Vector to multiply.
     * @return Resulting vector.
     * @throws std::invalid_argument if dimensions are incompatible.
     */
    VectorAlgebra<T> operator*(const VectorAlgebra<T> &vec) const;
};

template <typename T>
bool Matrix<T>::isValidate() const
{
    if (this->size() == 0)
        return true;

    size_t target = (*this)[0].size();
    for (const auto &row : *this)
    {
        if (row.size() != target)
            return false;
    }
    return true;
}

template <typename T>
Matrix<T> Matrix<T>::linearInterpolation(const Matrix<T> &m1, const Matrix<T> &m2, float ratio)
{
    if (m1.size() != m2.size())
        throw std::invalid_argument("Matrix row count mismatch");

    Matrix<T> result;
    result.reserve(m1.size());

    for (size_t i = 0; i < m1.size(); ++i)
        result.push_back(VectorAlgebra<T>::linearInterpolation(m1[i], m2[i], ratio));

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::projection(float fov_deg, float aspect_ratio, float near_plane, float far_plane)
{
    /*
      tan(θ) = opposite / adjacent
      so
      tan(fov / 2) = half_screen_height / near
      half_screen_height = near * tan(fov / 2)
      screen_height = 2 * near * tan(fov / 2)
      but instead to express a ratio between -1 and 1 we say
      scale = 2 / screen_height
      so
      scale = 2 / (2 * near * tan(fov / 2))
            = 1 / (near * tan(fov / 2))

      The role of X and Y scale is just to convert angles into [-1, 1] range — so the near value is not needed there anymore.
      y_scale = 1 / tan(fov / 2)
      x_scale = 1.0f / (tan(fov / 2) * aspect_ratio);

      for z:
      near * z_scale + z_translate = 0     → equation (1)
      far  * z_scale + z_translate = 1     → equation (2)
      Now subtract (1) from (2):
      (far - near) * z_scale = 1
      z_scale = 1 / (far - near)
      Plug into (1):
      near * (1 / (far - near)) + z_translate = 0
      z_translate = - near / (far - near)
      That’s what OpenGL used before
      But modern perspective projection includes a non-linear depth curve to make closer objects more precise in depth testing.
      The general form in OpenGL (for reversed Z and more precision) is often:
      z_scale     = far / (far - near)
      z_translate = - (far * near) / (far - near)
      This introduces nonlinearity, making:
      Depth buffer more precise near the camera
      Still gives 0 for near and 1 for far
    */
    float fov_rad = fov_deg * M_PI / 180.0f;       // Convert FOV to radians
    float tan_half_fov = std::tan(fov_rad / 2.0f); // tan(fov / 2)

    float x_scale = 1.0f / (tan_half_fov * aspect_ratio); // Horizontal FOV scaling (X axis). If X is far away, scale it less. If it's close, scale it more. Use FOV to compute the right factor.
    float y_scale = 1.0f / tan_half_fov;                  // Vertical FOV scaling (Y axis). If Y is far away, scale it less. If it's close, scale it more. Use FOV to compute the right factor.

    float z_range = far_plane - near_plane;                  // Depth range
    float z_scale = far_plane / z_range;                     // Z axis depth scale
    float z_translate = -(far_plane * near_plane) / z_range; // Shifts the z values to map the range [near, far] into [0, 1], while preserving depth precision

    /*
      We multiply a 4D vector {x, y, z, w} by this matrix.
      This matrix applies the perspective projection.
      After the multiplication, OpenGL will divide all x', y', z' by w' (called the perspective divide).
      Here's what the math looks like (AFTER transpose, i.e. real order used):
        x' = x * x_scale
        y' = y * y_scale
        z' = z * z_scale + w * 1
        w' = z * z_translate
      Then:
        x_ndc = x' / w'
        y_ndc = y' / w'
        z_ndc = z' / w'
      This is how 3D becomes screen coordinates!
    */
    Matrix<T> perspective = Matrix<T>::from({{x_scale, 0, 0, 0},
                                             {0, y_scale, 0, 0},
                                             {0, 0, z_scale, z_translate},
                                             {0, 0, 1, 0}});
    return perspective;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
    if (this->empty())
        return {};

    size_t rows = this->size();
    size_t cols = (*this)[0].size();

    Matrix<T> result;
    result.reserve(cols);

    for (size_t x = 0; x < cols; ++x)
    {
        VectorAlgebra<T> newRow;
        newRow.reserve(rows);
        for (size_t y = 0; y < rows; ++y)
            newRow.push_back((*this)[y][x]);
        result.push_back(newRow);
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::getIdentityMatrix() const
{
    Matrix<T> result;
    result.reserve(this->size());
    for (size_t y = 0; y < this->size(); y++)
    {
        result.emplace_back(this->size(), T(0));
        result[y][y] = T(1);
    }
    return result;
}

template <typename T>
unsigned int Matrix<T>::rank() const
{
    auto isNotZero = [](T t)
    { return t != T(0); };

    unsigned int rank = 0;
    Matrix<T> buffer = this->rowEchelon();

    for (const auto &row : buffer)
    {
        if (std::find_if(row.begin(), row.end(), isNotZero) != row.end())
            rank++;
    }

    return rank;
}

template <typename T>
Matrix<T> Matrix<T>::rowEchelon() const
{
    Matrix<T> buffer = *this;
    for (size_t y = 0; y < buffer.size(); y++)
    {
        for (size_t x = 0; x < buffer[y].size(); x++)
        {
            if (buffer[y][x] == T(0))
            {
                for (size_t yp = y + 1; yp < buffer.size(); yp++)
                {
                    if (buffer[yp][x] != T(0))
                    {
                        std::swap(buffer[yp], buffer[y]);
                        break;
                    }
                }
            }
            if (buffer[y][x] != T(0))
            {
                buffer[y] = buffer[y] / buffer[y][x];
                for (size_t yp = 0; yp < buffer.size(); yp++)
                {
                    if (yp != y)
                        buffer[yp] = buffer[yp] - (buffer[yp][x] * buffer[y]);
                }
                break;
            }
        }
    }
    return buffer;
}

template <typename T>
Matrix<T> Matrix<T>::rowEchelonNormalize(int &swapCount) const
{
    Matrix<T> buffer = *this;
    for (size_t y = 0; y < buffer.size(); y++)
    {
        for (size_t x = 0; x < buffer[y].size(); x++)
        {
            if (buffer[y][x] == T(0))
            {
                for (size_t yp = y + 1; yp < buffer.size(); yp++)
                {
                    if (buffer[yp][x] != T(0))
                    {
                        std::swap(buffer[yp], buffer[y]);
                        swapCount++;
                        break;
                    }
                }
            }
            if (buffer[y][x] != T(0))
            {
                for (size_t yp = 0; yp < buffer.size(); yp++)
                    if (yp != y)
                        buffer[yp] = buffer[yp] - ((buffer[yp][x] / buffer[y][x]) * buffer[y]);
                break;
            }
        }
    }
    return buffer;
}

template <typename T>
Matrix<T> Matrix<T>::invert() const
{
    if (!this->isSquare())
        throw std::logic_error("Cannot invert non-square matrix.");
    Matrix<T> identity = this->getIdentityMatrix();
    Matrix<T> result = *this;
    for (size_t y = 0; y < this->size(); y++)
        result[y].insert(result[y].end(), identity[y].begin(), identity[y].end());
    result = result.rowEchelon();
    for (size_t y = 0; y < this->size(); y++)
        result[y].assign(result[y].begin() + (result[y].size() / 2), result[y].end());
    return result;
}

template <typename T>
T Matrix<T>::determinant() const
{
    if (this->size() != (*this)[0].size())
        throw std::invalid_argument("Matrix must be square");

    int swapCount = 0;
    Matrix<T> upper = this->rowEchelonNormalize(swapCount);

    T det = 1;
    for (size_t i = 0; i < upper.size(); ++i)
        det *= upper[i][i];

    return (swapCount % 2 == 0) ? det : T(-1) * det;
}

template <typename T>
T Matrix<T>::trace()
{
    if (!isSquare())
        throw std::logic_error("Cannot take trace of this matrix. It is not square");
    if (this->empty())
        throw std::logic_error("Cannot take trace of an empty matrix.");
    T result = 0;
    size_t offset = 0;
    for (const auto &row : *this)
    {
        result += row[offset];
        offset++;
    }
    return result;
}

template <typename T>
bool Matrix<T>::isSquare() const
{
    size_t targetSize = this->size();
    for (const auto &row : *this)
    {
        if (row.size() != targetSize)
            return false;
    }
    return true;
}

template <typename T>
Matrix<T> Matrix<T>::from(const std::vector<VectorAlgebra<T>> &rows)
{
    Matrix<T> result;
    result.insert(result.end(), rows.begin(), rows.end());

    if (!result.isValidate())
        throw std::invalid_argument("All rows must be the same size to form a valid matrix");

    return result;
}

template <typename T>
template <typename Scalar>
Matrix<T> Matrix<T>::operator*(const Scalar scalar) const
{
    Matrix<T> result;
    for (auto row : (*this))
        result.push_back(row * scalar);
    return result;
}

template <typename Scalar, typename T>
Matrix<T> operator*(Scalar scalar, const Matrix<T> matrix)
{
    return matrix * scalar;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const
{
    if (this->empty() || other.empty() || (*this)[0].size() != other.size())
        throw std::invalid_argument("Incompatible matrix dimensions for multiplication");

    Matrix<T> otherBuffer = other.transpose();
    Matrix<T> result;
    result.reserve(this->size());

    for (const auto &rowA : *this)
    {
        VectorAlgebra<T> newRow;
        newRow.reserve(otherBuffer.size());

        for (const auto &rowB : otherBuffer)
            newRow.push_back(rowA.dotProduct(rowB));

        result.push_back(newRow);
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const
{
    if (this->empty() || (*this).size() != other.size() || (*this)[0].size() != other[0].size())
        throw std::invalid_argument("Matrix sizes mismatch");
    Matrix<T> result;
    for (size_t i = 0; i < (*this).size(); i++)
        result.push_back((*this)[i] + other[i]);
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const
{
    if (this->empty() || (*this).size() != other.size() || (*this)[0].size() != other[0].size())
        throw std::invalid_argument("Matrix sizes mismatch");
    Matrix<T> result;
    for (size_t i = 0; i < (*this).size(); i++)
        result.push_back((*this)[i] - other[i]);
    return result;
}

template <typename T>
VectorAlgebra<T> Matrix<T>::operator*(const VectorAlgebra<T> &vec) const
{
    if (this->empty() || (*this)[0].size() != vec.size())
        throw std::invalid_argument("Matrix and vector size mismatch");

    VectorAlgebra<T> result;
    result.reserve(this->size());

    for (const auto &row : *this)
        result.push_back(row.dotProduct(vec));

    return result;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    os << "[\n";
    for (const auto &row : matrix)
    {
        os << "  " << row << "\n";
    }
    os << "]";
    return os;
}
