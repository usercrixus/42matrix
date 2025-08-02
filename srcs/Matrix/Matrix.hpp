#pragma once

#include "VectorAlgebra.hpp"

template <typename T>
class Matrix : public std::vector<VectorAlgebra<T>>
{
private:
    bool validate();
    bool valid = false;
public:
    using std::vector<VectorAlgebra<T>>::vector;

    static Matrix<T> linearInterpolation(const Matrix<T> &m1, const Matrix<T> &m2, float ratio);

    Matrix<T> transpose() const;
    T trace();
    bool isSquare();

    Matrix<T> operator*(const Matrix<T> &other) const;
    VectorAlgebra<T> operator*(const VectorAlgebra<T> &vec) const;
};

template <typename T>
bool Matrix<T>::validate()
{
    if (this->size() == 0)
        return true;

    size_t target = (*this)[0].size();
    for (const auto &row : source)
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
bool Matrix<T>::isSquare()
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
