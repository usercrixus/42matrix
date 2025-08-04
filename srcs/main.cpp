#include "Matrix/VectorAlgebra.hpp"
#include "Matrix/Matrix.hpp"
#include "iostream"
#include <complex>

int main(int argc, char const *argv[])
{
	std::cout << "---- Exercice 00 ----" << std::endl;
	try
	{
		VectorAlgebra<int> vector1 = {0, 1, 2, 3, 4, 5};
		VectorAlgebra<int> vector2 = {0, 1, 2, 3, 4, 5};
		VectorAlgebra<int> vector3 = {0, 1, 2, 3, 4};
		std::cout << vector1 * vector2 << std::endl;
		std::cout << vector1 + vector2 << std::endl;
		std::cout << vector1 - vector2 << std::endl;
		std::cout << vector1 * 2 << std::endl;
		std::cout << vector1 - vector3 << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 01 ----" << std::endl;
		VectorAlgebra<float> e1 = {1, 0, 0};
		VectorAlgebra<float> e2 = {0, 1, 0};
		VectorAlgebra<float> e3 = {0, 0, 1};
		VectorAlgebra<float> v1 = {1, 2, 3};
		VectorAlgebra<float> v2 = {0, 10, -100};

		std::cout << VectorAlgebra<float>::linearCombinaison({e1, e2, e3}, {10, -2, 0.5}) << std::endl;
		std::cout << VectorAlgebra<float>::linearCombinaison({v1, v2}, {10, -2}) << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 02 ----" << std::endl;
		std::cout << VectorAlgebra<float>::linearInterpolation({0}, {1}, 0) << std::endl;
		std::cout << VectorAlgebra<float>::linearInterpolation({0}, {1}, 1) << std::endl;
		std::cout << VectorAlgebra<float>::linearInterpolation({0}, {1}, 0.5) << std::endl;
		std::cout << VectorAlgebra<float>::linearInterpolation({21}, {42}, 0.3) << std::endl;
		std::cout << VectorAlgebra<float>::linearInterpolation({2, 1}, {4, 2}, 0.3) << std::endl;
		std::cout << Matrix<float>::linearInterpolation(
						 Matrix<float>::from({{2, 1}, {3, 4}}),
						 Matrix<float>::from({{20, 10}, {30, 40}}), 0.5)
				  << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 03 ----" << std::endl;
		VectorAlgebra<float> v1 = {0, 0};
		VectorAlgebra<float> v2 = {1, 1};
		VectorAlgebra<float> v3 = {-1, 6};
		VectorAlgebra<float> v4 = {3, 2};
		std::cout << v1.dotProduct(v2) << std::endl;
		std::cout << v2.dotProduct(v2) << std::endl;
		std::cout << v3.dotProduct(v4) << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		/**
		 * can be used for normalization, or compare 2 vector difference, or stabilisation
		 */
		std::cout << "---- Exercice 04 ----" << std::endl;
		VectorAlgebra<float> v1 = {0, 0, 0};
		VectorAlgebra<float> v2 = {1, 2, 3};
		VectorAlgebra<float> v3 = {-1, -2};
		std::cout << v1.normManhattan() << " " << v1.normEuclidean() << " " << v1.normSupremum() << std::endl;
		std::cout << v2.normManhattan() << " " << v2.normEuclidean() << " " << v2.normSupremum() << std::endl;
		std::cout << v3.normManhattan() << " " << v3.normEuclidean() << " " << v3.normSupremum() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 05 ----" << std::endl;
		/**
		 * v1, v1	1	Same direction (0°)
		 * v1, v2	0	Perpendicular (90°)
		 * v3, v4	-1	Opposite direction (180°)
		 * v5, v6	1	Same direction (scaled)
		 * v7, v8	≈0.97	Small angle, similar direction
		 */
		VectorAlgebra<float> v1 = {1, 0};
		VectorAlgebra<float> v2 = {0, 1};
		VectorAlgebra<float> v3 = {-1, 1};
		VectorAlgebra<float> v4 = {1, -1};
		VectorAlgebra<float> v5 = {2, 1};
		VectorAlgebra<float> v6 = {4, 2};
		VectorAlgebra<float> v7 = {1, 2, 3};
		VectorAlgebra<float> v8 = {4, 5, 6};
		std::cout << VectorAlgebra<float>::angleCos(v1, v1) << std::endl;
		std::cout << VectorAlgebra<float>::angleCos(v1, v2) << std::endl;
		std::cout << VectorAlgebra<float>::angleCos(v3, v4) << std::endl;
		std::cout << VectorAlgebra<float>::angleCos(v5, v6) << std::endl;
		std::cout << VectorAlgebra<float>::angleCos(v7, v8) << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		/**
		 * usefull to find a perpendicular vector to 2 3d vector (plane)
		 */
		std::cout << "---- Exercice 06 ----" << std::endl;
		VectorAlgebra<float> v1 = {0, 0, 1};
		VectorAlgebra<float> v2 = {1, 0, 0};
		VectorAlgebra<float> v3 = {1, 2, 3};
		VectorAlgebra<float> v4 = {4, 5, 6};
		VectorAlgebra<float> v5 = {4, 2, -3};
		VectorAlgebra<float> v6 = {-2, -5, 16};
		std::cout << VectorAlgebra<float>::crossProduct(v1, v2) << std::endl;
		std::cout << VectorAlgebra<float>::crossProduct(v3, v4) << std::endl;
		std::cout << VectorAlgebra<float>::crossProduct(v5, v6) << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 07 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from({{1, 0}, {0, 1}});
		VectorAlgebra<float> v1 = {4, 2};
		Matrix<float> m2 = Matrix<float>::from({{2, 0}, {0, 2}});
		Matrix<float> m3 = Matrix<float>::from({{2, -2}, {-2, 2}});
		Matrix<float> m4 = Matrix<float>::from({{2, 1}, {4, 2}});
		Matrix<float> m5 = Matrix<float>::from({{3, -5}, {6, 8}});
		VectorAlgebra<float> v5 = {4, 2, -3};
		VectorAlgebra<float> v6 = {-2, -5, 16};
		std::cout << m1 * v1 << std::endl;
		std::cout << m2 * v1 << std::endl;
		std::cout << m3 * v1 << std::endl;
		std::cout << m1 * m1 << std::endl;
		std::cout << m1 * m4 << std::endl;
		std::cout << m5 * m4 << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 08 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from({{1, 0}, {0, 1}});
		Matrix<float> m2 = Matrix<float>::from({{2, -5, 0}, {4, 3, 7}, {-2, 3, 4}});
		Matrix<float> m3 = Matrix<float>::from({{-2, -8, 4}, {1, -23, 4}, {0, 6, 4}});
		std::cout << m1.trace() << std::endl;
		std::cout << m2.trace() << std::endl;
		std::cout << m3.trace() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 09 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from({{-2, -8, 4}, {1, -23, 4}, {0, 6, 4}});
		std::cout << m1 << std::endl;
		std::cout << m1.transpose() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 10 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from(
			{{1, 0, 0},
			 {0, 1, 0},
			 {0, 0, 1}});
		Matrix<float> m2 = Matrix<float>::from(
			{
				{1, 2},
				{3, 4},
			});
		Matrix<float> m3 = Matrix<float>::from(
			{
				{1, 2},
				{2, 4},
			});
		Matrix<float> m4 = Matrix<float>::from(
			{{8, 5, -2, 4, 28},
			 {4, 2.5, 20, 4, -4},
			 {8, 5, 1, 4, 17}});
		Matrix<float> m5 = Matrix<float>::from(
			{
				{0, 1},
				{1, 0},
			});
		std::cout << m1.rowEchelon() << std::endl;
		std::cout << m2.rowEchelon() << std::endl;
		std::cout << m3.rowEchelon() << std::endl;
		std::cout << m4.rowEchelon() << std::endl;
		std::cout << m5.rowEchelon() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 11 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from(
			{
				{1, -1},
				{-1, 1},
			});
		Matrix<float> m2 = Matrix<float>::from(
			{{2, 0, 0},
			 {0, 2, 0},
			 {0, 0, 2}});
		Matrix<float> m3 = Matrix<float>::from(
			{{8, 5, -2},
			 {4, 7, 20},
			 {7, 6, 1}});

		Matrix<float> m4 = Matrix<float>::from(
			{{8, 5, -2, 4},
			 {4, 2.5, 20, 4},
			 {8, 5, 1, 4},
			 {28, -4, 17, 1}});

		std::cout << m1.determinant() << std::endl;
		std::cout << m2.determinant() << std::endl;
		std::cout << m3.determinant() << std::endl;
		std::cout << m4.determinant() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 12 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from(
			{
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1},
			});
		Matrix<float> m2 = Matrix<float>::from(
			{
				{2, 0, 0},
				{0, 2, 0},
				{0, 0, 2},
			});
		Matrix<float> m3 = Matrix<float>::from(
			{
				{8, 5, -2},
				{4, 7, 20},
				{7, 6, 1},
			});

		std::cout << m1.invert() << std::endl;
		std::cout << m2.invert() << std::endl;
		std::cout << m3.invert() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 13 ----" << std::endl;
		Matrix<float> m1 = Matrix<float>::from(
			{
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1},
			});
		Matrix<float> m2 = Matrix<float>::from(
			{
				{1, 2, 0, 0},
				{2, 4, 0, 0},
				{-1, 2, 1, 1},
			});
		Matrix<float> m3 = Matrix<float>::from(
			{
				{8, 5, -2},
				{4, 7, 20},
				{7, 6, 1},
				{21, 18, 7},
			});

		std::cout << m1.rank() << std::endl;
		std::cout << m2.rank() << std::endl;
		std::cout << m3.rank() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 14 ----" << std::endl;
		float fov = 90.0f;		 // Field of view in degrees
		float ratio = 16.0f / 9; // Aspect ratio (e.g. 1920 / 1080)
		float near = 0.1f;		 // The closest distance you can see. 10cm
		float far = 100.0f;		 // The farthest distance you can see. 100 meters
		Matrix<float> m1 = Matrix<float>::projection(fov, ratio, near, far);
		std::cout << m1 << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	try
	{
		std::cout << "---- Exercice 15 (Complex numbers) ----" << std::endl;
		using Complex = std::complex<float>;

		// ----- Vector Tests -----
		VectorAlgebra<Complex> v1 = {{1.0f, 2.0f}, {3.0f, -4.0f}};
		VectorAlgebra<Complex> v2 = {{0.0f, -1.0f}, {1.0f, 1.0f}};

		std::cout << "v1: " << v1 << std::endl;
		std::cout << "v2: " << v2 << std::endl;
		std::cout << "v1 + v2: " << v1 + v2 << std::endl;
		std::cout << "v1 - v2: " << v1 - v2 << std::endl;
		std::cout << "v1 * v2 (element-wise): " << v1 * v2 << std::endl;
		std::cout << "v1 * scalar (2.0f + i): " << v1 * Complex(2.0f, 1.0f) << std::endl;
		std::cout << "dotProduct(v1, v2): " << v1.dotProduct(v2) << std::endl;
		std::cout << "normEuclidean(v1): " << v1.normEuclidean() << std::endl;
		std::cout << "normManhattan(v1): " << v1.normManhattan() << std::endl;
		std::cout << "normSupremum(v1): " << v1.normSupremum() << std::endl;
		std::cout << "v1 sum: " << v1.sum() << std::endl;

		auto interpolated = VectorAlgebra<Complex>::linearInterpolation(v1, v2, 0.5f);
		std::cout << "Interpolation v1/v2 @0.5: " << interpolated << std::endl;

		VectorAlgebra<Complex> coef = {{0.5f, 0.0f}, {0.5f, 0.0f}};
		auto combo = VectorAlgebra<Complex>::linearCombinaison({v1, v2}, coef);
		std::cout << "Linear combinaison: " << combo << std::endl;

		VectorAlgebra<Complex> v3 = {{1.0f, 0.0f}, {0.0f, 1.0f}, {1.0f, 1.0f}};
		VectorAlgebra<Complex> v4 = {{0.0f, 1.0f}, {1.0f, 0.0f}, {1.0f, -1.0f}};
		std::cout << "Cross product: " << VectorAlgebra<Complex>::crossProduct(v3, v4) << std::endl;

		// ----- Matrix Tests -----
		Matrix<Complex> m1 = Matrix<Complex>::from({{{1.0f, 2.0f}, {0.0f, 0.0f}},
													{{0.0f, 0.0f}, {1.0f, -2.0f}}});
		VectorAlgebra<Complex> vec = {{1.0f, 1.0f}, {1.0f, -1.0f}};

		std::cout << "\n-- Matrix operations --" << std::endl;
		std::cout << "Matrix m1:\n"
				  << m1 << std::endl;
		std::cout << "Matrix m1 * vec: " << m1 * vec << std::endl;
		std::cout << "Transpose of m1:\n"
				  << m1.transpose() << std::endl;
		std::cout << "Trace of m1: " << m1.trace() << std::endl;
		std::cout << "Determinant of m1: " << m1.determinant() << std::endl;
		std::cout << "Rank of m1: " << m1.rank() << std::endl;
		std::cout << "Inverted m1:\n"
				  << m1.invert() << std::endl;

		Matrix<Complex> identity = m1.getIdentityMatrix();
		std::cout << "Identity matrix:\n"
				  << identity << std::endl;
		std::cout << "m1 * identity:\n"
				  << m1 * identity << std::endl;

		auto interpolatedM = Matrix<Complex>::linearInterpolation(m1, identity, 0.5f);
		std::cout << "Interpolated matrix (m1 vs identity):\n"
				  << interpolatedM << std::endl;

		// ----- Projection -----
		std::cout << "\n-- Projection Matrix --" << std::endl;
		auto proj = Matrix<Complex>::projection(90.0f, 1.0f, 0.1f, 100.0f);
		std::cout << "Projection matrix (complex context):\n"
				  << proj << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
	}

	return 0;
}
