#ifndef __steer_steering_sim_params_h_
#define __steer_steering_sim_params_h_
#ifndef NO_STEER

class simulationParameters {

	public:

		double sim_pressure_min;
		double sim_pressure_max;
		double sim_velocity_min;
		double sim_velocity_max;
		double sim_stress_min;
		double sim_stress_max;
		int sim_time_step;
		double sim_time;
		int sim_cycle;
		int sim_n_inlets;
		double* sim_inlet_avg_vel;
		double sim_mouse_pressure;
		double sim_mouse_stress;

		XDR xdr_sim_params;
		char* sim_params;
		u_int sim_params_bytes;

		simulationParameters();
		~simulationParameters();
		char* pack();
		void collectGlobalVals();

};

#endif // NO_STEER

#endif//__steer_steering_sim_params_h_
