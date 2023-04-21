#include "CVRP_1.55_approximation.h"
#include <chrono>

// Constant pi
static const itv pi = kv::constants<itv>::pi();

// Define the empty interval as [inf,-inf]
// It works correctly with itv::hull function.
static const itv EmptyInterval(std::numeric_limits<itv::base_type>::infinity(), -std::numeric_limits<itv::base_type>::infinity());

// Template function to compute the integral of one (area) or distance (IntegralDistance) over a right triangle
// RightTriangle<Area>(x, y) = A0(x, y) in Definition 20
// RightTriangle<IntegralDistance>(x, y) = A1(x, y) in Definition 20
template<typename F> itv RightTriangle(const itv& x, const itv& y) {
	if constexpr (std::is_same_v<F, Area>) {
		return x * y / 2.;
	}
	else if constexpr (std::is_same_v<F, IntegralDistance>) {
		// If x can be close to zero, bound the IntegralDistance by (MaxDistance * Area) in order to avoid dividing by zero
		if (overlap(x, itv(-0.001,0.001))) {
			itv z(sqrt(abs(x * x + y * y)) * x * y / 2.);
			return itv::hull(z, -z);
		}
		// The formula in Lemma 21
		return x * (x * x * asinh(y / abs(x)) + y * sqrt(abs(x * x + y * y))) / 6.;
	}
}

// Function to compute the average distance between (x ,y) and a random point in [0, 1]^2
// AverageDistance(x, y) = g1(O), where O=(x, y)
itv AverageDistance(const itv& x, const itv& y) {
	// The formula in Theorem 22
	return RightTriangle<IntegralDistance>(x, y) + RightTriangle<IntegralDistance>(y, x) + RightTriangle<IntegralDistance>(y, 1. - x)
		+ RightTriangle<IntegralDistance>(1. - x, y) + RightTriangle<IntegralDistance>(1. - x, 1. - y)
		+ RightTriangle<IntegralDistance>(1. - y, 1. - x) + RightTriangle<IntegralDistance>(1. - y, x)
		+ RightTriangle<IntegralDistance>(x, 1. - y);
}

// Template function to compute the integral over a disk sector with unit radius and central angle theta
// DiskSector<Area>(theta) = theta / 2
// DiskSector<IntegralDistance>(theta) = theta / 3
template<typename F> itv DiskSector(const itv& theta) {
	if constexpr (std::is_same_v<F, Area>)
		return theta / 2.;
	else if constexpr (std::is_same_v<F, IntegralDistance>)
		return theta / 3.;
}

// Template function to compute the integral over a disk segment with unit radius and apothem h
// DiskSegment<Area>(h) = B0(h) in Definition 23
// DiskSegment<IntegralDistance>(h) = B1(h) in Definition 23
template <typename F> itv DiskSegment(const itv& h) {
	itv a = min(itv(1.), max(itv(-1.), h));
	itv b = sqrt(abs(1. - a * a));
	itv theta = 2. * (pi - acos(a));
	// The formula in Lemma 24
	return DiskSector<F>(theta) + 2. * RightTriangle<F>(a, b);
}

// Template function to compute the integral over the intersection of a unit disk and two half-planes with parameters h1, h2
// DiskDoubleSegment<Area>(h1, h2) = C0(h1, h2) in Definition 25
// DiskDoubleSegment<IntegralDistance>(h1, h2) = C1(h1, h2) in Definition 25
template <typename F> itv DiskDoubleSegment(const itv& h1, const itv& h2) {
	itv result(EmptyInterval);
	// The formula in Lemma 26
	// Check if (h1, h2) may be outside the unit disk
	if ((h1 * h1 + h2 * h2).upper() >= 1.) {
		// Check if the first case in Theorem 26 may happen
		if (h1.lower() <= 0. && h2.lower() <= 0.) {
			result = itv::hull(result, itv(0.));
		}
		// Check if the second case in Theorem 26 may happen
		if (h1.upper() >= 0. && h2.lower() <= 0.) {
			result = itv::hull(result, DiskSegment<F>(h2));
		}
		// Check if the third case in Theorem 26 may happen
		if (h1.lower() <= 0. && h2.upper() >= 0.) {
			result = itv::hull(result, DiskSegment<F>(h1));
		}
		// Check if the fourth case in Theorem 26 may happen
		if (h1.upper() >= 0. && h2.upper() >= 0.) {
			result = itv::hull(result, DiskSegment<F>(h1) + DiskSegment<F>(h2) - DiskSector<F>(2. * pi));
		}
	}
	// Check if (h1, h2) may be inside the unit disk
	if (overlap(h1 * h1 + h2 * h2, itv(0., 1.))) {
		itv a1 = min(itv(1.), max(itv(-1.), h1));
		itv a2 = min(itv(1.), max(itv(-1.), h2));
		itv b1 = sqrt(abs(1. - a1 * a1));
		itv b2 = sqrt(abs(1. - a2 * a2));
		itv theta = pi / 2. + asin(a1) + asin(a2);
		// If so, we also need to consider the fifth case in Theorem 26
		result = itv::hull(result, DiskSector<F>(theta) + RightTriangle<F>(a1, b1) + RightTriangle<F>(a2, b2) + RightTriangle<F>(a1, a2) + RightTriangle<F>(a2, a1));
	}
	// Let the result be the convex hull of all possible cases
	return result;
}

// Template function to compute the functions D0 and D1 in Definition 27
// OverlapRegion<Area>(x, y, R) = D0(x, y, R) in Definition 27
// OverlapRegion<IntegralDistance>(x, y, R) = D1(x, y, R) in Definition 27
template <typename F> itv OverlapRegion(const itv& x, const itv& y, const itv& R) {
	// The formula in Theorem 28
	return DiskDoubleSegment<F>((1. - x) / R, (1. - y) / R) - DiskDoubleSegment<F>((1. - x) / R, -y / R)
		- DiskDoubleSegment<F>(-x / R, (1. - y) / R) + DiskDoubleSegment<F>(-x / R, -y / R);
}

int main() {
	// An itv variable x is a closed inteval [x.lower(), x.upper()] represented by a pair of C++ double variables
	// Operations on itv variables are precision-guaranteed 
	// We want interval result1 to cover all possible values of g2(O') - 31/48 * g1(O') 
	//  and interval result2 to cover all possible values of g3(O') - 31/48 
	itv result1(EmptyInterval);
	itv result2(EmptyInterval);
	itv a, b;
	itv R, g1_val, g2_val, g3_val, D0_val;
	const itv c1(itv(31.) / 48.);
	const itv c2(itv(1.) / 500.);
	// Timing start
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i <= 2371; i++) {
		for (int j = i; j <= 2371; j++) {
			// Verify the inequalities for every O' = (a, b) in the epsilon net N
			// We don't worry about 0.5 and 0.75 here because they have finite binary representations:
			//  they have no rounding errors as C++ double
			a = 0.5 + c2 * i;
			b = 0.5 + c2 * j;
			// We have g1_val = g1(O'), g2_val = g2(O') and g3_val = g3(O')
			g1_val = AverageDistance(a, b);
			// The formulas in Theorem 29
			R = g1_val * 0.75;
			D0_val = OverlapRegion<Area>(a, b, R);
			g2_val = R - R * R * R * D0_val + R * R * R * OverlapRegion<IntegralDistance>(a, b, R);
			g3_val = 1 - R * R * D0_val;
			result1 = itv::hull(result1, g2_val - g1_val * c1);
			result2 = itv::hull(result2, g3_val - c1);
		}
	}
	// Takes about half an hour on my laptop
	auto end = std::chrono::high_resolution_clock::now();
	std::cout.precision(16);
	std::cout << std::boolalpha;
	std::cout << result1.lower() << "," << result2.lower() << std::endl;
	std::cout << "g2(O') - 31/48 * g1(O') >= 0.0025 is " << (result1 >= "0.0025") << std::endl;
	std::cout << "g3(O') - 31/48 >= 0.0096 is " << (result2 >= "0.0096") << std::endl;
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Running time：" << elapsed.count() << " seconds" << std::endl; 
	return 0;
}
