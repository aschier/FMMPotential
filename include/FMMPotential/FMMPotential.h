#pragma once

#include <vector>
#include <valarray>
#include <memory>
#include <Eigen/Eigen>

namespace FMMPotential {
	typedef unsigned int uint;
	struct Triangle {
		size_t p1, p2, p3;
		size_t operator[](int i) const {
			if(i == 0) {
				return p1;
			} else if(i == 1) {
				return p2;
			} else if(i==2) {
				return p3;
			}
			assert(false);
			return 0;
		}
	};

	typedef std::vector<std::valarray<double>> VVertices;
	typedef std::vector<Triangle> VTriangles;

	class FMMPotential {
	private:
		class impl;
		std::unique_ptr<impl> pImpl;
	public:
		FMMPotential(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit);
		FMMPotential(const VVertices &vertices, const VTriangles &triangles, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
		FMMPotential(const std::vector<std::vector<std::valarray<double>>> &triangle_points, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
		~FMMPotential();
		void evaluateField(Eigen::MatrixX3d  &field, const Eigen::MatrixX3d  &particles);
		void evaluatePotentials(Eigen::VectorXd &potentials, const Eigen::MatrixX3d  &particles);
		void evaluatePotentialsAndField(Eigen::VectorXd &potentials, Eigen::MatrixX3d  &field, const Eigen::MatrixX3d  &particles);
	};
}
