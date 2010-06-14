#ifndef __steering_sim_params_h
#define __steering_sim_params_h

#ifndef NO_STEER

// Class for encapsulating the parameters of the simulation for sending to the client
class simulationParameters {

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

  // TODO: After refactoring, these methods should be made protected (when things are in namespaces).
  public:
    // Setter methods
    void set_Min_Sim_Pressure(double new_min_pressure);
    void set_Max_Sim_Pressure(double new_max_pressure);
    void set_Min_Sim_Velocity(double new_min_velocity);
    void set_Max_Sim_Velocity(double new_max_velocity);
    void set_Min_Sim_Stress(double new_min_stress);
    void set_Max_Sim_Stress(double new_max_stress);

    void set_Sim_Inlets(int n_inlets);

    void set_Sim_Mouse_Pressure(double new_mouse_pressure);
    void set_Sim_Mouse_Stress(double new_mouse_stress);    

  public:
    simulationParameters();
    ~simulationParameters();

    // Return a char array of the packed parameters
    char* pack();
    // TODO: Remove this method.
    void collectGlobalVals();
    // Number of bytes required to pack the simulation parameters
    u_int getPackedSizeInBytes();

    // Accessor methods
    double get_Min_Sim_Pressure();
    double get_Max_Sim_Pressure();
    double get_Min_Sim_Velocity();
    double get_Max_Sim_Velocity();
    double get_Min_Sim_Stress();
    double get_Max_Sim_Stress();
};

#endif // NO_STEER

#endif // __steering_sim_params_h
