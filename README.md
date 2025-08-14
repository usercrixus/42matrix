# Matrix & VectorAlgebra — Modern C++ Linear Algebra (Header-only)

A small, header-only linear algebra toolkit featuring:

- **`VectorAlgebra<T>`** — a `std::vector<T>` with element-wise ops, dot/cross products, norms, interpolation, and linear combinations.
- **`Matrix<T>`** — a matrix built on `std::vector<VectorAlgebra<T>>` with transpose, row-echelon, rank, determinant, inverse, multiplication, trace, interpolation, and a perspective projection helper.

Works with real types (`int`, `float`, `double`) and complex (`std::complex<float>`, `std::complex<double>`). Complex support follows standard linear algebra conventions (Hermitian dot, real-valued norms).

---

## Table of contents
- [Install & Build](#install--build)
- [Hello world](#hello-world)
- [Design notes](#design-notes)
- [VectorAlgebra API](#vectoralgebra-api)
- [Matrix API](#matrix-api)
- [Complex numbers](#complex-numbers)
- [Errors & exceptions](#errors--exceptions)
- [FAQ & tips](#faq--tips)
- [License](#license)

---

## Install & Build

This library is **header-only**. Just add the two headers to your project:

```
Matrix/
 ├─ VectorAlgebra.hpp
 └─ Matrix.hpp
```

Include them:

```cpp
#include "Matrix/VectorAlgebra.hpp"
#include "Matrix/Matrix.hpp"
```

### Build example (g++)
```bash
g++ -std=c++20 -Wall -Wextra -Werror main.cpp -o demo
./demo
```

---

## Hello world

```cpp
#include "Matrix/VectorAlgebra.hpp"
#include "Matrix/Matrix.hpp"
#include <iostream>

int main() {
    VectorAlgebra<float> a = {1, 2, 3};
    VectorAlgebra<float> b = {4, 5, 6};

    std::cout << "a + b = " << (a + b) << "\n";
    std::cout << "a · b = " << a.dotProduct(b) << "\n";
    std::cout << "||a||2 = " << a.normEuclidean() << "\n";

    Matrix<float> M = Matrix<float>::from({{1, 2}, {3, 4}});
    VectorAlgebra<float> x = {1, 1};
    std::cout << "M * x = " << (M * x) << "\n";
}
```

Example output:
```
a + b = [5, 7, 9]
a · b = 32
||a||2 = 3.74166
M * x = [3, 7]
```

---

## Design notes

- `VectorAlgebra<T>` **inherits publicly** from `std::vector<T>` to keep standard vector behaviors (iterators, indexing, etc.) and add math methods/operators.
- `Matrix<T>` **inherits privately** from `std::vector<VectorAlgebra<T>>` to control the public surface (exposes `operator[]`, `begin`, `end`).
- All operations do **dimension checks** and throw exceptions on mismatch.

---

## VectorAlgebra API

### Type & construction
```cpp
template <typename T>
class VectorAlgebra : public std::vector<T> {
public:
    using std::vector<T>::vector; // inherit constructors
    using ScalarType = decltype(std::abs(T{})); // real scalar type for norms
};
```

Create like a normal vector:
```cpp
VectorAlgebra<float> v = {1, 2, 3};
VectorAlgebra<int>   u(5, 0);   // [0,0,0,0,0]
```

### Static utilities
```cpp
static VectorAlgebra<T> linearCombinaison(const std::vector<VectorAlgebra<T>>& vects, const VectorAlgebra<T>& coef);
static VectorAlgebra<T> linearInterpolation(const VectorAlgebra<T>& a, const VectorAlgebra<T>& b, float ratio);
static float angleCos(const VectorAlgebra<T>& a, const VectorAlgebra<T>& b);
```

### Vector metrics
```cpp
ScalarType normManhattan() const; // ∑ |x_i|
ScalarType normEuclidean() const; // sqrt(∑ |x_i|^2)
ScalarType normSupremum() const;  // max_i |x_i|
```

### Reductions & products
```cpp
T sum();
T dotProduct(const VectorAlgebra<T>& v) const;
```

### Element-wise operators
```cpp
VectorAlgebra<T> operator+(const VectorAlgebra<T>&) const;
VectorAlgebra<T> operator-(const VectorAlgebra<T>&) const;
VectorAlgebra<T> operator*(const VectorAlgebra<T>&) const;
VectorAlgebra<T> operator/(const VectorAlgebra<T>&) const;
template <typename Scalar>
VectorAlgebra<T> operator*(Scalar s) const;
template <typename Scalar>
VectorAlgebra<T> operator/(Scalar s) const;
template <typename Scalar, typename U>
VectorAlgebra<U> operator*(Scalar s, const VectorAlgebra<U>& v);
template <typename Scalar, typename U>
VectorAlgebra<U> operator/(Scalar s, const VectorAlgebra<U>& v);
```

### Cross product (3D only)
```cpp
static VectorAlgebra<T> crossProduct(const VectorAlgebra<T>& a, const VectorAlgebra<T>& b);
```

---

## Matrix API

### Construction
```cpp
static Matrix<T> from(const std::vector<VectorAlgebra<T>>& rows);
bool isSquare() const;
```

### Basic ops
```cpp
Matrix<T> transpose() const;
VectorAlgebra<T> operator*(const VectorAlgebra<T>& x) const;
Matrix<T> operator*(const Matrix<T>& B) const;
```

### Row-echelon & rank
```cpp
Matrix<T> rowEchelon() const;
unsigned int rank() const;
Matrix<T> rowEchelonNormalize(int& swapCount) const;
```

### Determinant, trace, inverse
```cpp
T determinant() const;
T trace();
Matrix<T> invert() const;
```

### Interpolation
```cpp
static Matrix<T> linearInterpolation(const Matrix<T>& A, const Matrix<T>& B, float ratio);
```

### Perspective projection helper
```cpp
static Matrix<T> projection(float fov_deg, float aspect_ratio, float near_plane, float far_plane);
```

---

## Complex numbers

- **Dot product:** `a.dotProduct(b)` computes `∑ conj(a_i) * b_i`.
- **Norms:** return real scalar types (`ScalarType`) via `std::abs`/`std::norm`.
- **Angle cosine:** intentionally **throws** for complex vectors.
- All matrix routines work with `std::complex`.

---

## Errors & exceptions

Most precondition failures throw `std::invalid_argument`.

`std::logic_error` is thrown for conceptual errors.

---

## FAQ & tips

- **Why inherit from `std::vector`?** Convenience & zero-overhead reuse.
- **Performance:** Not optimized for heavy workloads; prefer Eigen/Blaze for that.
- **Numerical stability:** Only basic pivoting implemented.
- **Printing:** `operator<<` defined for both vectors and matrices.

---

## License

MIT License