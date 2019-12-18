#pragma once
#include <float.h>
#include <cassert>

typedef unsigned int uint;

struct rule {
	const int num_points;
	const double(*points)[2];
	const double *weights;
};

const static double rule2_points[3][2] = {
	{1./3, 1./3},
	{1./3, 1./6},
	{1./6, 1./3}
};
const static double rule2_weights[3] = {1./3, 1./3, 1./3};
const static rule rule2{3, rule2_points, rule2_weights};

const static double rule4_points[6][2]{
	{0.44594849091597, 0.44594849091597},
	{0.44594849091597, 0.10810301816807},
	{0.10810301816807, 0.44594849091597},
	{0.09157621350977, 0.09157621350977},
	{0.09157621350977, 0.81684757298046},
	{0.81684757298046, 0.09157621350977}
};
const static double rule4_weights[6]{
	0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532
};
const static rule rule4{6, rule4_points, rule4_weights};

const static double rule6_points[12][2]{
	{0.24928674517091, 0.24928674517091},
	{0.24928674517091, 0.50142650965818},
	{0.50142650965818, 0.24928674517091},
	{0.06308901449150, 0.06308901449150},
	{0.06308901449150, 0.87382197101700},
	{0.87382197101700, 0.06308901449150},
	{0.31035245103378, 0.63650249912140},
	{0.63650249912140, 0.05314504984482},
	{0.05314504984482, 0.31035245103378},
	{0.63650249912140, 0.31035245103378},
	{0.31035245103378, 0.05314504984482},
	{0.05314504984482, 0.63650249912140}
};
const static double rule6_weights[12]{
	0.11678627572638, 0.11678627572638, 0.11678627572638,
	0.05084490637021, 0.05084490637021, 0.05084490637021,
	0.08285107561837, 0.08285107561837, 0.08285107561837,
	0.08285107561837, 0.08285107561837, 0.08285107561837
};
const static rule rule6{12, rule6_points, rule6_weights};
const static rule rules[3]{rule2, rule4, rule6};

#define DEFAULT_QUADRATURE_DEGREE 6

inline rule get_quad_rule(uint polynomial_degree) {
	short selected_rule = 0;
	if(polynomial_degree < 3) {
		selected_rule = 0;
	} else if(polynomial_degree < 5) {
		selected_rule = 1;
	} else if(polynomial_degree < 7) {
		selected_rule = 2;
	} else {
		assert(false);
	}
	return rules[selected_rule];
}

/**
* Quadratic Gauss quadrature rule for triangles. Gauss points are at 1/6 1/6, 1/6 2/3 and 2/3 1/6.
* @param a The first triangle vertex as array with 3 coordinates.
* @param b The first triangle vertex as array with 3 coordinates.
* @param c The first triangle vertex as array with 3 coordinates.
* @param function The function to integrate, which maps a array
*        with 3 values to the function value at the point.
* @returns The approximated integral divided by the triangle area.
*/
template<typename Scalar>
inline Scalar gaussTriangleQuad(const Scalar a[3], const Scalar b[3], const Scalar c[3], unsigned short quad_degree, Scalar function(Scalar[3], void *data), void *data = nullptr) {
	Scalar result = 0.0;
	rule quad_rule = get_quad_rule(quad_degree);
	Scalar point[3];
	for(short i = 0; i < quad_rule.num_points; i++) {
		for(short j = 0; j < 3; j++) {
			point[j] = quad_rule.points[i][0] * a[j] + quad_rule.points[i][1] * b[j] + (1.0 - (quad_rule.points[i][0] + quad_rule.points[i][1])) * c[j];
		}
		result += quad_rule.weights[i] * function(point, data);
	}

	// 2*area is the absolute value of the jacobian |(d(x,y)/(d(xi, eta))| of the nodal shape function transformation.
	// We return the result divided by area (by not multiplying with it), as we do not have the area available here and
	// do not want to calculate it here.
	return result * 2;
}

template<typename Scalar>
inline void gaussTriangleQuad(Scalar result[3], const Scalar a[3], const Scalar b[3], const Scalar c[3], unsigned short quad_degree, void function(Scalar[3], Scalar[3], void *data), void *data = nullptr) {
	result[0] = 0.0; result[1] = 0.0; result[2] = 0.0;
	rule quad_rule = get_quad_rule(quad_degree);
	Scalar point[3];
	Scalar function_result[3];
	for(short i = 0; i < quad_rule.num_points; i++) {
		for(short j = 0; j < 3; j++) {
			point[j] = quad_rule.points[i][0] * a[j] + quad_rule.points[i][1] * b[j] + (1 - (quad_rule.points[i][0] + quad_rule.points[i][1])) * c[j];
		}
		function(function_result, point, data);
		for(short j = 0; j < 3; j++) {
			// 2*area is the absolute value of the jacobian |(d(x,y)/(d(xi, eta))| of the nodal shape function transformation.
			// We return the result divided by area (by not multiplying with it), as we do not have the area available here and
			// do not want to calculate it here.
			result[j] += 2 * quad_rule.weights[i] * function_result[j];
		}
	}
}
