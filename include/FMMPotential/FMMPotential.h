#pragma once

#include <vector>
#include <array>
#include <memory>
#include <Eigen/Eigen>

#ifdef WITH_MISHMESH
#include <MishMesh/TriMesh.h>
#endif

namespace FMMPotential {
	typedef unsigned int uint;
	typedef std::array<size_t, 3> Triangle;

	class FMMPotential {
	private:
		class impl;
		std::unique_ptr<impl> pImpl;
	public:
		FMMPotential(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit);
		FMMPotential(const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
		FMMPotential(const std::vector<Eigen::Matrix3d> &triangle_points, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
#ifdef WITH_MISHMESH
		FMMPotential(const MishMesh::TriMesh &mesh, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
#endif
		~FMMPotential();
		void evaluateField(Eigen::MatrixX3d  &field, const Eigen::MatrixX3d  &particles);
		void evaluatePotentials(Eigen::VectorXd &potentials, const Eigen::MatrixX3d  &particles);
		void evaluatePotentialsAndField(Eigen::VectorXd &potentials, Eigen::MatrixX3d  &field, const Eigen::MatrixX3d  &particles);
	};
}
