#ifndef NO_STEER

// Class for encapsulating the parameters of the simulation for sending to the client
class simulationParameters {

  // TODO Make this class the way simulation parameters are stored internally in the simulation - otherwise this is a bit useless
  private:
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
    double sim_mouse_pressure;
    double sim_mouse_stress;

    XDR xdr_sim_params;
    char* sim_params;
    
  public:
    simulationParameters();
    ~simulationParameters();

    // Return a char array of the packed parameters
    char* pack();
    void collectGlobalVals();
    // Number of bytes required to pack the simulation parameters
    u_int getPackedSizeInBytes();
};

#endif // NO_STEER
