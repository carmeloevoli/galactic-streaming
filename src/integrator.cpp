#include "integrator.hpp"

using namespace integrators;

double integrator_base::step_adaptive(utils::Grid<double> &waves,  utils::Grid<double> &particles,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_waves,
        std::function<utils::Grid<double>(utils::Grid<double> &, utils::Grid<double> &)> &rhs_particles,
        double time, double & del_t, double epsilon0) {

	// Store previous step
	utils::Grid<double> waves_old = waves;
	utils::Grid<double> particles_old = particles;

	bool step_accepted = false;

	int iter=0;
	do {
		// For an error estimate, do full timestep AND two half timesteps

		// One full step
		double time_full = time;
		time_full = step(waves, particles, rhs_waves, rhs_particles, time_full, del_t);
		// store results
		utils::Grid<double> waves_full= waves;
		utils::Grid<double> particles_full = particles;

		// reset fields
		waves = waves_old;
		particles = particles_old;

		// Two half steps
		double time_half = time;
		time_half = step(waves, particles, rhs_waves, rhs_particles, time_half, 0.5*del_t);
		time_half = step(waves, particles, rhs_waves, rhs_particles, time_half, 0.5*del_t);

		// Compare results
		utils::Grid<double> diff_waves = waves_full - waves;
		double eps_waves = diff_waves.get_RMS()/(pow(2.,scheme_order) - 1.);
		utils::Grid<double> diff_particles = particles_full - particles;
		double eps_particles = diff_particles.get_RMS()/(pow(2.,scheme_order) - 1.);

		double epsilon = std::max(eps_waves, eps_particles);

		iter++;

		if(epsilon < epsilon0) {
			step_accepted = true;

			// Compute new delta t (possibly larger than before)
			del_t = beta*del_t*pow(epsilon0/epsilon, 1./(scheme_order + 1.));
			time = time_half;

		} else {
			// if step not accepted reset to previous step and repeat
			waves = waves_old;
			particles = particles_old;
			// Compute new delta t (smaller thann before)
			del_t = beta*del_t*pow(epsilon0/epsilon, 1./(scheme_order));
		}


	} while (!step_accepted);

	if(m_verbosity>3) {
		std::cout << " Result after " << iter << " steps with " << del_t << "\n";
	}

	return time;

}
