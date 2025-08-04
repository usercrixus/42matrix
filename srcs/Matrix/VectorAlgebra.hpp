#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <complex>

template <typename T>
class VectorAlgebra : public std::vector<T>
{
private:
public:
	// Use base constructor
	using std::vector<T>::vector;

	using ScalarType = decltype(std::abs(T{}));
	static VectorAlgebra<T> linearCombinaison(const std::vector<VectorAlgebra<T>> &vects, const VectorAlgebra<T> &coef);
	static VectorAlgebra<T> linearInterpolation(const VectorAlgebra<T> &vect1, const VectorAlgebra<T> &vect2, float ratio);
	static float angleCos(const VectorAlgebra<T> &vect1, const VectorAlgebra<T> &vect2);
	static VectorAlgebra<T> crossProduct(const VectorAlgebra<T> &vect1, const VectorAlgebra<T> &vect2);
	ScalarType normManhattan() const;
	ScalarType normEuclidean() const;
	ScalarType normSupremum() const;

	T sum();
	T dotProduct(const VectorAlgebra<T> &other) const;

	VectorAlgebra<T> operator+(const VectorAlgebra<T> &other) const;
	VectorAlgebra<T> operator-(const VectorAlgebra<T> &other) const;
	VectorAlgebra<T> operator*(const VectorAlgebra<T> &other) const;
	template <typename Scalar>
	VectorAlgebra<T> operator*(const Scalar) const;
	VectorAlgebra<T> operator/(const VectorAlgebra<T> &other) const;
	template <typename Scalar>
	VectorAlgebra<T> operator/(Scalar value) const;
};

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::linearCombinaison(const std::vector<VectorAlgebra<T>> &vects, const VectorAlgebra<T> &coef)
{
	if (vects.size() != coef.size())
		throw std::invalid_argument("Vector size mismatch");

	size_t size = vects[0].size();
	for (auto &vec : vects)
	{
		if (vec.size() != size)
			throw std::invalid_argument("All vectors must be the same size");
	}

	VectorAlgebra<T> result(size, T(0));
	for (size_t i = 0; i < vects.size(); i++)
		result = result + (vects[i] * coef[i]);

	return result;
}

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::linearInterpolation(const VectorAlgebra<T> &vect1, const VectorAlgebra<T> &vect2, float ratio)
{
	if (vect1.size() != vect2.size())
		throw std::invalid_argument("vect1 and vect2 should have the same size");
	if (ratio < 0 || ratio > 1)
		throw std::invalid_argument("ratio should be in [0, 1]");

	return ((1 - ratio) * vect1) + (ratio * vect2);
}

template <typename T>
float VectorAlgebra<T>::angleCos(const VectorAlgebra<T> &vect1, const VectorAlgebra<T> &vect2)
{
	if constexpr (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value)
		throw std::logic_error("angleCos is undefined for complex vectors.");
	T dot = vect1.dotProduct(vect2);
	T normA = vect1.normEuclidean();
	T normB = vect2.normEuclidean();

	if (normA == 0 || normB == 0)
		throw std::invalid_argument("Cannot compute cosine with zero-length vector");

	return dot / (normA * normB);
}

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::crossProduct(const VectorAlgebra<T> &vect1, const VectorAlgebra<T> &vect2)
{
	if (vect1.size() != 3 || vect2.size() != 3)
		throw std::invalid_argument("Cross product is only defined for 3D vectors");

	return {vect1[1] * vect2[2] - vect1[2] * vect2[1], vect1[2] * vect2[0] - vect1[0] * vect2[2], vect1[0] * vect2[1] - vect1[1] * vect2[0]};
}

template <typename T>
typename VectorAlgebra<T>::ScalarType VectorAlgebra<T>::normManhattan() const
{
	typename VectorAlgebra<T>::ScalarType sum = 0;
	for (const auto &val : *this)
		sum += std::abs(val);
	return sum;
}

template <typename T>
typename VectorAlgebra<T>::ScalarType VectorAlgebra<T>::normEuclidean() const
{
	typename VectorAlgebra<T>::ScalarType sum = 0;
	for (const auto &val : *this)
		sum += std::norm(val);
	return std::sqrt(sum);
}

template <typename T>
typename VectorAlgebra<T>::ScalarType VectorAlgebra<T>::normSupremum() const
{
	typename VectorAlgebra<T>::ScalarType maxAbs = std::abs((*this)[0]);
	for (const auto &val : *this)
	{
		typename VectorAlgebra<T>::ScalarType  magnitude = std::abs(val);
		if (magnitude > maxAbs)
			maxAbs = magnitude;
	}
	return maxAbs;
}

template <typename T>
T VectorAlgebra<T>::sum()
{
	T result = 0;
	for (auto &value : (*this))
		result += value;
	return result;
}

template <typename T>
T VectorAlgebra<T>::dotProduct(const VectorAlgebra<T> &other) const
{
	if (this->size() != other.size())
		throw std::invalid_argument("Dot product: vector size mismatch");

	T result = 0;
	for (size_t i = 0; i < this->size(); ++i)
	{
		if constexpr (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value)
			result += std::conj((*this)[i]) * other[i];
		else
			result += (*this)[i] * other[i];
	}

	return result;
}

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::operator+(const VectorAlgebra<T> &other) const
{
	if (this->size() != other.size())
		throw std::invalid_argument("Vector size mismatch");

	VectorAlgebra<T> result;
	result.reserve(this->size());

	for (size_t i = 0; i < this->size(); ++i)
		result.push_back((*this)[i] + other[i]);

	return result;
}

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::operator-(const VectorAlgebra<T> &other) const
{
	if (this->size() != other.size())
		throw std::invalid_argument("Vector size mismatch");

	VectorAlgebra<T> result;
	result.reserve(this->size());

	for (size_t i = 0; i < this->size(); ++i)
		result.push_back((*this)[i] - other[i]);

	return result;
}

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::operator*(const VectorAlgebra<T> &other) const
{
	if (this->size() != other.size())
		throw std::invalid_argument("Vector size mismatch");

	VectorAlgebra<T> result;
	result.reserve(this->size());

	for (size_t i = 0; i < this->size(); ++i)
		result.push_back((*this)[i] * other[i]);

	return result;
}

template <typename T>
template <typename Scalar>
VectorAlgebra<T> VectorAlgebra<T>::operator*(Scalar value) const
{
	VectorAlgebra<T> result;
	result.reserve(this->size());

	for (size_t i = 0; i < this->size(); ++i)
		result.push_back((*this)[i] * value);

	return result;
}

template <typename Scalar, typename T>
VectorAlgebra<T> operator*(Scalar value, const VectorAlgebra<T> &vec)
{
	return vec * value;
}

template <typename T>
VectorAlgebra<T> VectorAlgebra<T>::operator/(const VectorAlgebra<T> &other) const
{
	if (this->size() != other.size())
		throw std::invalid_argument("Vector size mismatch");

	VectorAlgebra<T> result;
	result.reserve(this->size());

	for (size_t i = 0; i < this->size(); ++i)
	{
		if (other[i] == 0)
			throw std::invalid_argument("cant divide by zero");
		result.push_back((*this)[i] / other[i]);
	}

	return result;
}

template <typename T>
template <typename Scalar>
VectorAlgebra<T> VectorAlgebra<T>::operator/(Scalar value) const
{
	if (value == Scalar(0))
		throw std::invalid_argument("cant divide by zero");
	VectorAlgebra<T> result;
	result.reserve(this->size());

	for (size_t i = 0; i < this->size(); ++i)
		result.push_back((*this)[i] / value);

	return result;
}

template <typename Scalar, typename T>
VectorAlgebra<T> operator/(Scalar value, const VectorAlgebra<T> &vec)
{
	VectorAlgebra<T> result;
	result.reserve(vec.size());

	for (size_t i = 0; i < vec.size(); ++i)
	{
		if (vec[i] == 0)
			throw std::invalid_argument("cant divide by zero");
		result.push_back((value / vec[i]));
	}

	return result;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const VectorAlgebra<T> &vec)
{
	os << "[";
	for (size_t i = 0; i < vec.size(); ++i)
	{
		os << vec[i];
		if (i != vec.size() - 1)
			os << ", ";
	}
	os << "]";
	return os;
}