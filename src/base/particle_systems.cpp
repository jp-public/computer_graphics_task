#include "particle_systems.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <chrono>

using namespace std;
using namespace FW;

namespace {

	inline Vec3f fGravity(float mass) {
		return Vec3f(0, -9.8f * mass, 0);
	}

	// force acting on particle at pos1 due to attached to pos2 at the other end
	inline Vec3f fSpring(const Vec3f& pos1, const Vec3f& pos2, float k, float rest_length) {
		// YOUR CODE HERE (R2)
		Vec3f d = pos1 - pos2;
		float d_len = d.length();
		Vec3f fspring = -k * (d_len - rest_length) * (d / d_len);
		return fspring;
	}

	inline Vec3f fDrag(const Vec3f& v, float k) {
		// YOUR CODE HERE (R2)
		return -k * v;
	}

	inline Vec3f fWind(boolean wind_) {
		Vec3f wind;
		if (wind_) {
			srand(time(NULL));
			float r1 = ((float)rand() / (float)(RAND_MAX)) * 2 -1.f;
			float r2 = ((float)rand() / (float)(RAND_MAX)) * 2 -1.f;
			float r3 = ((float)rand() / (float)(RAND_MAX)) * 2 - 1.f;
			wind = (r1, r2, r3);
		}
		return wind;
	}

} // namespace

void SimpleSystem::reset() {
	current_state_ = State(1, Vec3f(0, radius_, 0));
}

State SimpleSystem::evalF(const State& state) const {
	State f(1, Vec3f(-state[0].y, state[0].x, 0));
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H
// using the implicit Euler method, the simple system should converge towards origin -- as opposed to the explicit Euler, which diverges outwards from the origin.
void SimpleSystem::evalJ(const State&, SparseMatrix& result, bool initial) const {
	if (initial) {
		result.coeffRef(1, 0) = 1.0f;
		result.coeffRef(0, 1) = -1.0f;
	}
}
#endif

Points SimpleSystem::getPoints() {
	return Points(1, current_state_[0]);
}

Lines SimpleSystem::getLines() {
	static const auto n_lines = 50u;
	auto l = Lines(n_lines * 2);
	const auto angle_incr = 2 * FW_PI / n_lines;
	for (auto i = 0u; i < n_lines; ++i) {
		l[2 * i] = l[2 * i + 1] =
			Vec3f(radius_ * FW::sin(angle_incr * i), radius_ * FW::cos(angle_incr * i), 0);
	}
	rotate(l.begin(), l.begin() + 1, l.end());
	return l;
}

void SpringSystem::reset() {
	const auto start_pos = Vec3f(0.1f, -0.5f, 0.0f);
	const auto spring_k = 30.0f;
	const auto rest_length = 0.5f;
	current_state_ = State(4);
	// YOUR CODE HERE (R2)
	// Set the initial state for a particle system with one particle fixed
	// at origin and another particle hanging off the first one with a spring.
	// Place the second particle initially at start_pos.
	current_state_[0] = Vec3f(0.f, 0.f, 0.f);
	current_state_[1] = Vec3f(0.f, 0.f, 0.f);
	current_state_[2] = start_pos;
	current_state_[3] = Vec3f(0.f, 0.f, 0.f);
	spring_ = Spring(0, 1, spring_k, rest_length);
}

State SpringSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	State f(4);
	// YOUR CODE HERE (R2)
	// Return a derivative for the system as if it was in state "state".
	// You can use the fGravity, fDrag and fSpring helper functions for the forces.
	Vec3f g = fGravity(mass);
	Vec3f drag = fDrag(state[3], drag_k);
	auto s = spring_;
	Vec3f springf = fSpring(state[2], state[0], s.k, s.rlen);
	f[2] = state[3];
	f[3] = (g+drag+springf)/mass;
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

// This is a very useful read for the Jacobians of the spring forces. It deals with spring damping as well, we don't do that -- our drag is simply a linear damping of velocity (that results in some constants in the Jacobian).
// http://blog.mmacklin.com/2012/05/04/implicitsprings/

void SpringSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.5f;
	const auto mass = 1.0f;
	// EXTRA: Evaluate the Jacobian into the 'result' matrix here. Only the free end of the spring should have any nonzero values related to it.
}
#endif

Points SpringSystem::getPoints() {
	auto p = Points(2);
	p[0] = current_state_[0]; p[1] = current_state_[2];
	return p;
}

Lines SpringSystem::getLines() {
	auto l = Lines(2);
	l[0] = current_state_[0]; l[1] = current_state_[2];
	return l;
}

void PendulumSystem::reset() {
	const auto spring_k = 1000.0f;
	const auto start_point = Vec3f(0);
	const auto end_point = Vec3f(0.05, -1.5, 0);
	current_state_ = State(2 * n_);
	springs_.clear();
	
	// YOUR CODE HERE (R4)
	// Set the initial state for a pendulum system with n_ particles
	// connected with springs into a chain from start_point to end_point with uniform intervals.
	// The rest length of each spring is its length in this initial configuration.
	float x_step = 0.05f / (n_-1);
	float y_step = -1.5f / (n_-1);
	for (auto i = 0u; i < 2*n_; ++i) {
		if (i % 2 == 0) {
			current_state_[i] = Vec3f(0.f + (i/2) * x_step, 0.f + (i/2) * y_step, 0.f);
		} else {
			current_state_[i] = Vec3f(0.f, 0.f, 0.f);		
		}
	}
	float rest_length = end_point.length() / (n_-1);
	for (auto i = 0u; i < n_-1; i++) {
		springs_.push_back(Spring(i, i+1, spring_k, rest_length));
	}
}

State PendulumSystem::evalF(const State& state) const {
	const auto drag_k = 0.5f;
	const auto mass = 0.5f;
	auto f = State(2 * n_);
	// YOUR CODE HERE (R4)
	// As in R2, return a derivative of the system state "state".
	Vec3f springf;
	for (const auto& s : springs_) {
		springf = fSpring(state[s.i1 * 2], state[s.i2 * 2], s.k, s.rlen);
		f[s.i1 * 2 + 1] += springf;
		springf = fSpring(state[s.i2 * 2], state[s.i1 * 2], s.k, s.rlen);
		f[s.i2 * 2 + 1] += springf;
	}
	Vec3f g = fGravity(mass);
	Vec3f drag;
	for (auto i = 2u; i < 2 * n_; ++i) {
		if (i % 2 == 0) {
			f[i] = state[i + 1];
		}
		else {
			drag = fDrag(state[i], drag_k);
			f[i] = (f[i] + g + drag) / mass;
		}
	}
	f[1] = 0;
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void PendulumSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {

	const auto drag_k = 0.5f;
	const auto mass = 0.5f;

	// EXTRA: Evaluate the Jacobian here. Each spring has an effect on four blocks of the matrix -- both of the positions of the endpoints will have an effect on both of the velocities of the endpoints.
}
#endif


Points PendulumSystem::getPoints() {
	auto p = Points(n_);
	for (auto i = 0u; i < n_; ++i) {
		p[i] = current_state_[i * 2];
	}
	return p;
}

Lines PendulumSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(current_state_[2 * s.i1]);
		l.push_back(current_state_[2 * s.i2]);
	}
	return l;
}

void ClothSystem::reset() {
	const auto spring_k = 300.0f;
	const auto width = 1.5f, height = 1.5f; // width and height of the whole grid
	current_state_ = State(2 * x_*y_);
	// YOUR CODE HERE (R5)
	// Construct a particle system with a x_ * y_ grid of particles,
	// connected with a variety of springs as described in the handout:
	// structural springs, shear springs and flex springs.
	springs_.clear();
	wind_ = false;
	unsigned int n = x_ * y_;
	float step = 1.5f / (x_ - 1);
	for (auto i = 0u; i < x_; ++i) {
		for (auto j = 0u; j < y_; ++j) {
			current_state_[2*i * x_ + 2*j] = Vec3f((i * step)-0.75f, 0.f, j * -step);
		}
	}
	float rest_length = 1.5f / (x_-1);
	for (auto i = 0u; i < x_; i++) {
		for (auto j = 0u; j < y_; j++) {
			if (j < y_ - 1) {
				springs_.push_back(Spring(i * x_ + j, i * x_ + j + 1, spring_k, rest_length));
			}
			if (i < x_ - 1) {
				springs_.push_back(Spring(i * x_ + j, i * x_ + j + x_, spring_k, rest_length));
			}
			if (j < y_ - 2) {
				springs_.push_back(Spring(i * x_ + j, i * x_ + j + 2, spring_k, rest_length*2));
			}
			if (i < x_ - 2) {
				springs_.push_back(Spring(i * x_ + j, i * x_ + j + 2*x_, spring_k, rest_length*2));
			}
			if (j < y_ - 1 && i < x_ - 1) {
				springs_.push_back(Spring(i * x_ + j, i * x_ + j + x_+1, spring_k, rest_length*sqrt(2)));
			}
			if (j > 0 && i < x_ - 1) {
				springs_.push_back(Spring(i * x_ + j, i * x_ + j + x_-1, spring_k, rest_length*sqrt(2)));
			}
		}
	}
}

void ClothSystem::removeSprings(const State& state) {
	for (int i = springs_.size() - 1; i >= 0; i--) {
		Spring s = springs_[i];
		float current_len = (state[s.i1 * 2] - state[s.i2 * 2]).length();
		float tearing_factor = 1.63f;
		float tearing_threshold = s.rlen * tearing_factor;
		if (current_len > tearing_threshold) {
			springs_.erase(springs_.begin() + i);	
		}

	}
}

State ClothSystem::evalF(const State& state) const {
	const auto drag_k = 0.08f;
	const auto n = x_ * y_;
	static const auto mass = 0.025f;
	auto f = State(2 * n);
	// YOUR CODE HERE (R5)
	// This will be much like in R2 and R4.
	Vec3f springf;
	for (const auto& s : springs_) {
		springf = fSpring(state[s.i1 * 2], state[s.i2 * 2], s.k, s.rlen);
		f[s.i1 * 2 + 1] += springf;
		springf = fSpring(state[s.i2 * 2], state[s.i1 * 2], s.k, s.rlen);
		f[s.i2 * 2 + 1] += springf;
	}
	Vec3f g = fGravity(mass);
	Vec3f drag;
	Vec3f wind; 
	
	for (auto i = 0u; i < 2 * n; ++i) {
		if (i % 2 == 0) {
			f[i] = state[i + 1];
		}
		else {
			drag = fDrag(state[i], drag_k);
			wind = fWind(wind_);
			f[i] = (f[i] + g + drag+wind) / mass;
		}
	}
	f[1] = 0;
	f[2*n - y_*2] = 0;
	return f;
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void ClothSystem::evalJ(const State& state, SparseMatrix& result, bool initial) const {
	const auto drag_k = 0.08f;
	static const auto mass = 0.025f;

	// EXTRA: Evaluate the Jacobian here. The code is more or less the same as for the pendulum.
}

#endif

Points ClothSystem::getPoints() {
	auto n = x_ * y_;
	auto p = Points(n);
	for (auto i = 0u; i < n; ++i) {
		p[i] = current_state_[2 * i];
	}
	return p;
}

Lines ClothSystem::getLines() {
	auto l = Lines();
	for (const auto& s : springs_) {
		l.push_back(current_state_[2 * s.i1]);
		l.push_back(current_state_[2 * s.i2]);
	}
	return l;
}
State FluidSystem::evalF(const State&) const {
	return State();
}

