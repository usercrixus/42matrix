#include "Matrix/VectorAlgebra.hpp"
#include "Matrix/Matrix.hpp"
#include "iostream"

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
		std::cout << Matrix<float>::linearInterpolation({{2, 1}, {3, 4}}, {{20, 10}, {30, 40}}, 0.5) << std::endl;
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
		Matrix<float> m1 = {{1, 0}, {0, 1}};
		VectorAlgebra<float> v1 = {4, 2};
		Matrix<float> m2 = {{2, 0}, {0, 2}};
		Matrix<float> m3 = {{2, -2}, {-2, 2}};
		Matrix<float> m4 = {{2, 1}, {4, 2}};
		Matrix<float> m5 = {{3, -5}, {6, 8}};
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
		Matrix<float> m1 = {{1, 0}, {0, 1}};
		Matrix<float> m2 = {{2, -5, 0}, {4, 3, 7}, {-2, 3, 4}};
		Matrix<float> m3 = {{-2, -8, 4}, {1, -23, 4}, {0, 6, 4}};
		std::cout << m1.trace() << std::endl;
		std::cout << m2.trace() << std::endl;
		std::cout << m3.trace() << std::endl;
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	return 0;
}
