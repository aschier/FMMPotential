#pragma once

#include <vector>
#include <Eigen/Eigen>

#ifndef EXAFMM_INDEXED_BODIES
#define EXAFMM_INDEXED_BODIES
#endif
#include "exafmm-3d/exafmm.h"

#include <FMMPotential/FMMPotential.h>

using namespace std;
typedef unsigned int uint;

namespace FMMPotential {
	class FMMPotential::impl {
	public:
		impl(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta, uint P, uint ncrit);
		impl(const VVertices &vertices, const VTriangles &triangles, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
		impl(const vector<vector<valarray<double>>> &triangle_points, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
#ifdef WITH_MISHMESH
		impl(const MishMesh::TriMesh &mesh, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64, uint quad_degree = 6);
#endif
		void evaluateField(Eigen::MatrixX3d  &field, const Eigen::MatrixX3d  &particles);
		void evaluatePotentials(Eigen::VectorXd &potentials, const Eigen::MatrixX3d  &particles);
		void evaluatePotentialsAndField(Eigen::VectorXd &potentials, Eigen::MatrixX3d  &field, const Eigen::MatrixX3d  &particles);

	private:
		void init(const Eigen::MatrixX3d &charge_points, const Eigen::VectorXd &charges, double theta = 0.4, uint P = 5, uint ncrit = 64);
		void evaluate(const Eigen::MatrixX3d &particles);

		// instance data
		exafmm::Bodies charge_bodies;
		exafmm::Cells charge_cells;
		// stored as cache to avoid unneeded reallocation
		exafmm::Bodies particle_bodies;
	};

	// Conversion functions

	/**
	 * Get a vector with the force at a body.
	 * @param body The body.
	 * @returns A vector with the force.
	 */
	inline Eigen::Vector3d body2force(const exafmm::Body &body) {
		return {body.F[0], body.F[1], body.F[2]};
	}

	/**
	 * Get the coordinate vector from a body as vector.
	 * @param body The body.
	 * @returns A vector with the coordinates.
	 */
	inline Eigen::Vector3d body2point(const exafmm::Body &body) {
		return {body.X[0], body.X[1], body.X[2]};
	}

	/**
	 * Get a vector with the potentials at a Bodies at list.
	 * @param bodies A Bodies list of the bodies to get the potential from.
	 * @returns a vector with one entry per body.
	 */
	inline Eigen::VectorXd bodies2potentials(const exafmm::Bodies &bodies) {
		Eigen::VectorXd potentials(bodies.size());
		for(size_t i = 0; i < bodies.size(); i++) {
			potentials[i] = bodies[i].p;
		}
		return potentials;
	}

	/**
	* Get a vector with the charges at a Bodies at list.
	* @param bodies A Bodies list of the bodies to get the charge from.
	* @returns a vector with one entry per body.
	*/
	inline Eigen::VectorXd bodies2charges(const exafmm::Bodies &bodies) {
		Eigen::VectorXd charges(bodies.size());
		for(size_t i = 0; i < bodies.size(); i++) {
			charges[i] = bodies[i].q;
		}
		return charges;
	}

	/**
	* Get the force vectors from a list of bodies
	* @param bodies A Bodies list of the bodies to get the forces from.
	* @returns a Nx3 matrix with the force vectors of the bodies
	*/
	inline Eigen::MatrixX3d bodies2field(const exafmm::Bodies &bodies) {
		Eigen::MatrixX3d field(bodies.size(), 3);
		for(size_t i = 0; i < bodies.size(); i++) {
			field(i, 0) = bodies[i].F[0];
			field(i, 1) = bodies[i].F[1];
			field(i, 2) = bodies[i].F[2];
		}
		return field;
	}

	/**
	* Convert a Nx3 matrix of particles given by their position to bodies with 0 charge.
	* @param[out] A Bodies list, with bodies initialized with the coordinates and 0 charge/force/potential.
	*/
	inline void particles2bodies(exafmm::Bodies &bodies, const Eigen::MatrixX3d &particles) {
		bodies.resize(particles.rows());
		for(int i = 0; i < particles.rows(); i++) {
			bodies[i].X[0] = particles(i, 0);
			bodies[i].X[1] = particles(i, 1);
			bodies[i].X[2] = particles(i, 2);
			bodies[i].q = 0.0;
			bodies[i].p = 0.0;
			bodies[i].F = 0.0;
			bodies[i].index = i;
		}
	}

	/**
	 * Convert a Nx3 matrix with particles given by their position to bodies with 0 charge.
	 * @param particles The particles.
	 * @returns A Bodies list, with bodies initialized with the coordinates and 0 charge/force/potential.
	 */
	inline exafmm::Bodies particles2bodies(const Eigen::MatrixX3d particles) {
		exafmm::Bodies bodies(particles.rows());
		particles2bodies(bodies, particles);
		return bodies;
	}
}
