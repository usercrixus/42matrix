#pragma once

#include "VectorAlgebra.hpp"

template <typename T>
class Matrix : private std::vector<VectorAlgebra<T>>
{
private:
    bool isValidate() const;
    bool isSquare() const;

public:
    using std::vector<VectorAlgebra<T>>::operator[];
    using std::vector<VectorAlgebra<T>>::begin;
    using std::vector<VectorAlgebra<T>>::end;

    static Matrix<T> from(const std::vector<VectorAlgebra<T>> &rows);
    static Matrix<T> linearInterpolation(const Matrix<T> &m1, const Matrix<T> &m2, float ratio);

    Matrix<T> transpose() const;
    Matrix<T> rowEchelon() const;
    Matrix<T> rowEchelonNormalize(int &swapCount) const;
    Matrix<T> invert() const;
    Matrix<T> getIdentityMatrix() const;
    unsigned int rank() const;
    T determinant() const;
    T trace();

    Matrix<T> operator*(const Matrix<T> &other) const;
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
        result[y][y] = 1;
    }
    return result;
}

template <typename T>
unsigned int Matrix<T>::rank() const
{
    auto isNotZero = [](T t) { return t != T(0); };

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
            if (buffer[y][x] == 0)
            {
                for (size_t yp = y + 1; yp < buffer.size(); yp++)
                {
                    if (buffer[yp][x] != 0)
                    {
                        std::swap(buffer[yp], buffer[y]);
                        break;
                    }
                }
            }
            if (buffer[y][x] != 0)
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
            if (buffer[y][x] == 0)
            {
                for (size_t yp = y + 1; yp < buffer.size(); yp++)
                {
                    if (buffer[yp][x] != 0)
                    {
                        std::swap(buffer[yp], buffer[y]);
                        swapCount++;
                        break;
                    }
                }
            }
            if (buffer[y][x] != 0)
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

    return (swapCount % 2 == 0) ? det : -1 * det;
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
