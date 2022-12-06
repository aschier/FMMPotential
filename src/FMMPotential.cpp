#include "FMMPotential/FMMPotential.h"

#include <vector>
#include <valarray>
#include <cassert>
#include <algorithm>
#include <iostream>

#ifndef EXAFMM_INDEXED_BODIES
#define EXAFMM_INDEXED_BODIES
#endif
#include "exafmm-3d/exafmm.h"
#include "exafmm-3d/build_tree.h"
#include "exafmm-3d/kernel.h"
#include "exafmm-3d/traverse_eager.h"

#include "quadrature.hpp"
#include "FMMPotential_impl.h"

#ifdef WITH_MISHMESH
#include <MishMesh/macros.h>
#include <MishMesh/utils.h>
#endif

#ifndef _OPENMP
Breaks on purpose. Your openmp is not working.
#endif

using namespace std;

namespace FMMPotential {
	/**
	 * Create a FMMPotential object from a vector of charge points.
	 * The initial potentials and forces on the given bodies will be calculated and used for evaluation on other particles.
	 * @param vertices A vector with the coordinates of the charge points.
	 * @param charges The charges of the points.
	 * @param theta Fast Multipole Theta parameter.
	 * @param P polynomial degree in the Fast Multipole Method.
	 * @param ncrit Size of the tree leafs in the Fast Multipole Method.
	 */
	FMMPotential::impl::impl(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit) {
		assert(charge_points.rows() == charges.rows());
		init(charge_points, charges, theta, P, ncrit);
	}
	/**
	 * Create a FMMPotential object from vectors of vertices, triangles and a valarray with the charges.
	 * The initial potentials and forces on the given bodies will be calculated and used for evaluation on other particles.
	 * @param vertices A vector of triangle vertices.
	 * @param triangles A vector of triangles used for integration.
	 * @param charges The total charges of the triangles.
	 * @param theta Fast Multipole Theta parameter.
	 * @param P polynomial degree in the Fast Multipole Method.
	 * @param ncrit Size of the tree leafs in the Fast Multipole Method.
	 * @param quad_degree Polynomial degree the quadrature rule is exact for. Possible values: 2, 4, 6.
	 *        Values 0, 1, 3, 5 will use the next higher degree of exactness and values bigger than 6 are
	 *        not supported.
	 */
	FMMPotential::impl::impl(const VVertices &vertices, const VTriangles &triangles, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit, uint quad_degree) {
		rule quad_rule = get_quad_rule(quad_degree);

		Eigen::MatrixX3d charge_points(quad_rule.num_points * triangles.size(), 3);
		Eigen::VectorXd point_charges(quad_rule.num_points * triangles.size());
		// Split the charge of a triangle into point charges at gauss quadrature points.
		// TODO: maybe this should be abstracted in some function in quadrature.h
		for(size_t t_idx = 0; t_idx < triangles.size(); t_idx++) {
			const valarray<double> &a = vertices[triangles[t_idx][0]];
			const valarray<double> &b = vertices[triangles[t_idx][1]];
			const valarray<double> &c = vertices[triangles[t_idx][2]];
			for(short i = 0; i < quad_rule.num_points; i++) {
				for(short j = 0; j < 3; j++) {
					charge_points(t_idx * quad_rule.num_points + i, j) = quad_rule.points[i][0] * a[j] + quad_rule.points[i][1] * b[j] + (1.0 - (quad_rule.points[i][0] + quad_rule.points[i][1])) * c[j];
				}
				// 2 because of the transformation between reference triangle and actual triangle, see quadrature.h
				point_charges[t_idx * quad_rule.num_points + i] = charges[t_idx] * quad_rule.weights[i] * 2;
			}
		}
		init(charge_points, point_charges, theta, P, ncrit);
	}

	/**
	 * Create a FMMPotential object from a vector of triangles given by their vertices and a valarray with the charges.
	 * The initial potentials and forces on the given bodies will be calculated and used for evaluation on other particles.
	 * @param triangle_points A vector of triangles given by their vertices.
	 * @param charges The total charges of the triangles.
	 * @param theta Fast Multipole Theta parameter.
	 * @param P polynomial degree in the Fast Multipole Method.
	 * @param ncrit Size of the tree leafs in the Fast Multipole Method.
	 * @param quad_degree Polynomial degree the quadrature rule is exact for. Possible values: 2, 4, 6.
	 *        Values 0, 1, 3, 5 will use the next higher degree of exactness and values bigger than 6 are
	 *        not supported.
	 */
	FMMPotential::impl::impl(const vector<vector<valarray<double>>> &triangle_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit, uint quad_degree) {
		rule quad_rule = get_quad_rule(quad_degree);

		Eigen::MatrixX3d charge_points(quad_rule.num_points * triangle_points.size(), 3);
		Eigen::VectorXd point_charges(quad_rule.num_points * triangle_points.size());
		// Split the charge of a triangle into point charges at gauss quadrature points.
		// TODO: maybe this should be abstracted in some function in quadrature.h
		for(size_t t_idx = 0; t_idx < triangle_points.size(); t_idx++) {
			const valarray<double> &a = triangle_points[t_idx][0];
			const valarray<double> &b = triangle_points[t_idx][1];
			const valarray<double> &c = triangle_points[t_idx][2];
			for(short i = 0; i < quad_rule.num_points; i++) {
				for(short j = 0; j < 3; j++) {
					charge_points(t_idx * quad_rule.num_points + i, j) = quad_rule.points[i][0] * a[j] + quad_rule.points[i][1] * b[j] + (1.0 - (quad_rule.points[i][0] + quad_rule.points[i][1])) * c[j];
				}
				point_charges[t_idx * quad_rule.num_points + i] = charges[t_idx] * quad_rule.weights[i];
			}
		}
		init(charge_points, point_charges, theta, P, ncrit);
	}

#ifdef WITH_MISHMESH
	FMMPotential::impl::impl(const MishMesh::TriMesh & mesh, const Eigen::VectorXd & charges, double theta, uint P, uint ncrit, uint quad_degree) {
		rule quad_rule = get_quad_rule(quad_degree);

		Eigen::MatrixX3d charge_points(quad_rule.num_points * mesh.n_faces(), 3);
		Eigen::VectorXd point_charges(quad_rule.num_points * mesh.n_faces());
		// Split the charge of a triangle into point charges at gauss quadrature points.
		for(int t_idx = 0; t_idx < mesh.n_faces(); t_idx++) {
			auto fh = mesh.face_handle(t_idx);
			auto points = MishMesh::face_points(mesh, fh);
			for(short i = 0; i < quad_rule.num_points; i++) {
				for(short j = 0; j < 3; j++) {
					charge_points(t_idx * quad_rule.num_points + i, j) = quad_rule.points[i][0] * points[0][j] + quad_rule.points[i][1] * points[1][j] + (1.0 - (quad_rule.points[i][0] + quad_rule.points[i][1])) * points[2][j];
				}
				// 2 because of the transformation between reference triangle and actual triangle, see quadrature.h
				point_charges[t_idx * quad_rule.num_points + i] = charges[t_idx] * quad_rule.weights[i] * 2;
			}
		}
		init(charge_points, point_charges, theta, P, ncrit);
	}
#endif

	/**
	 * Evaluate the potential of the charged bodies in the instance on the given particles.
	 * @param[in] particles A Nx3 matrix storing the particle positions in its rows.
	 * @note the result is stored in the instance variable particle_bodies and is not used directly.
	 */
	void FMMPotential::impl::evaluate(const Eigen::MatrixX3d &particles) {
		particles2bodies(particle_bodies, particles);
		exafmm::Cells particle_cells = exafmm::buildTree(particle_bodies);
		for(size_t i = 0; i < particle_cells.size(); i++) {
			particle_cells[i].L.resize(exafmm::NTERM);
			particle_cells[i].M.resize(exafmm::NTERM);
		}
		exafmm::horizontalPass(particle_cells, charge_cells);
		exafmm::downwardPass(particle_cells);
		std::sort(particle_bodies.begin(), particle_bodies.end()); // exafmm reorders the bodies
	}

	/**
	 * Get the field for a vector of particles given by their coordinates.
	 * @param[out] field The resulting field.
	 * @param[in] particles A Nx3 matrix storing the particle positions in its rows.
	 * @note If you need potentials and field, use evaluatePotentialsAndField to avoid double calculation.
	 */
	void FMMPotential::impl::evaluateField(Eigen::MatrixX3d &field, const Eigen::MatrixX3d &particles) {
		evaluate(particles);
		for(size_t i = 0; i < particle_bodies.size(); i++) {
			field(i, 0) = particle_bodies[i].F[0];
			field(i, 1) = particle_bodies[i].F[1];
			field(i, 2) = particle_bodies[i].F[2];
		}
	}

	/**
	 * Get the potentials for a vector of particles given by their coordinates.
	 * @param[out] potentials The resulting potentials.
	 * @param[in] particles The particles vector.
	 * @note If you need potentials and field, use evaluatePotentialsAndField to avoid double calculation.
	 */
	void FMMPotential::impl::evaluatePotentials(Eigen::VectorXd &potentials, const Eigen::MatrixX3d &particles) {
		evaluate(particles);
		for(uint i = 0; i < particle_bodies.size(); i++) {
			potentials[i] = particle_bodies[i].p;
		}
	}

	/**
	 * Get the potentials and field for a vector of particles given by their coordinates.
	 * @param[out] potentials The resulting potentials.
	 * @param[out] field The resulting potentials.
	 * @param[in] particles A Nx3 matrix storing the particle positions in its rows.
	 * @note The method is more efficient than calling evaluatePotentials and evaluateField separately,
	 *       as the bodies only need to be calculated once.
	 */
	void FMMPotential::impl::evaluatePotentialsAndField(Eigen::VectorXd &potentials, Eigen::MatrixX3d &field, const Eigen::MatrixX3d &particles) {
		evaluate(particles);
		for(uint i = 0; i < particle_bodies.size(); i++) {
			potentials[i] = particle_bodies[i].p;
		}
		for(size_t i = 0; i < particle_bodies.size(); i++) {
			field(i, 0) = particle_bodies[i].F[0];
			field(i, 1) = particle_bodies[i].F[1];
			field(i, 2) = particle_bodies[i].F[2];
		}
	}

	void FMMPotential::impl::init(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit) {
		if(ncrit < 8) {
			ncrit = 8;
			cerr << "NCRIT too small for FMM, increasing to 8" << endl;
		}
		exafmm::THETA = theta;
		exafmm::NCRIT = ncrit;
		exafmm::P = P;

		assert(charge_points.rows() == charges.rows());
		charge_bodies = exafmm::Bodies(charge_points.rows());
		for(int i = 0; i < charge_points.rows(); i++) {
			exafmm::vec3 X;
			for(short j = 0; j < 3; j++) { X[j] = charge_points(i, j); };
			charge_bodies[i] = exafmm::Body{X, charges[i], 0, 0, i};
		}

		exafmm::initKernel();
		charge_cells = exafmm::buildTree(charge_bodies);
		exafmm::upwardPass(charge_cells);
	};
	FMMPotential::FMMPotential(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit){
		pImpl = std::make_unique<FMMPotential::impl>(charge_points, charges, theta, P, ncrit);
	}
	FMMPotential::FMMPotential(const VVertices &vertices, const VTriangles &triangles, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit, uint quad_degree) {
		pImpl = std::make_unique<FMMPotential::impl>(vertices, triangles, charges, theta, P, ncrit, quad_degree);
	}
	FMMPotential::FMMPotential(const std::vector<std::vector<std::valarray<double>>>&triangle_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit, uint quad_degree) {
		pImpl = std::make_unique<FMMPotential::impl>(triangle_points, charges, theta, P, ncrit, quad_degree);
	}
#ifdef WITH_MISHMESH
	FMMPotential::FMMPotential(const MishMesh::TriMesh & mesh, const Eigen::VectorXd & charges, double theta, uint P, uint ncrit, uint quad_degree) {
		pImpl = std::make_unique<FMMPotential::impl>(mesh, charges, theta, P, ncrit, quad_degree);
	}
#endif
	FMMPotential::~FMMPotential() = default;

	void FMMPotential::evaluateField(Eigen::MatrixX3d &field, const Eigen::MatrixX3d &particles) {
		pImpl->evaluateField(field, particles);
	}
	void FMMPotential::evaluatePotentials(Eigen::VectorXd &potentials, const Eigen::MatrixX3d &particles) {
		pImpl->evaluatePotentials(potentials, particles);
	}
	void FMMPotential::evaluatePotentialsAndField(Eigen::VectorXd &potentials, Eigen::MatrixX3d &field, const Eigen::MatrixX3d &particles) {
		pImpl->evaluatePotentialsAndField(potentials, field, particles);
	}
}
