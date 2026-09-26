#pragma once
#include "types.hpp"
#include "constants.hpp"
#include "sim/inc/properties.hpp"
#include "sim/inc/ins.hpp"
#include "sim/inc/data_export.hpp"
#include "fc/inc/fc_sim_connector.hpp"
#include <array>
#include <memory>
#include <string>
#include <vector>
#include "renderer/earth_surface.hpp"

struct RocketStartState {
    Vec3 origin_r_eci;
    Quat origin_q_eci; // origin attitude
    Vec3 target_r_ecef;
};

// snapshot of rocket state at any given moment
struct RocketState {
    bool detonation_active = false;

    double t = 0;
    double mass = 0, fuel = 0;
    double length = 0, cm_dist = 0, engine_dist = 0, radius = 0;   // dims from nose
    double nose_length = 0;
    bool has_engine = true;

    Vec3 r{}, v{}, a{}, w{};
    Quat q_rocket{1, 0, 0, 0};
    Quat q_engine{1, 0, 0, 0};
    RocketStartState init{};
};

class Rocket {
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // public                                                                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    public:

    // counter to track once the rocket is dead how long it should stay existing before deleting itself
    double life_countdown = 30.0; // stays alive for n (sim) seconds before disappearing

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // setup                                                                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    Rocket(const std::string& name, double origin_latitude, double origin_longitude, double target_latitude,
           double target_longitude, const RocketProps& props, bool track_data, double export_interval);
    ~Rocket();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // getters                                                                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    RocketState get_state() const;
    const std::string& get_name() const { return name; }
    bool is_detonated() { return detonated; }
    int active_stage_idx() const { return active_idx; } // index the fc's stage array with this

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // setters                                                                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    void set_pos(const Vec3& pos) { r = pos; } // set absolute position
    void set_orientation(const Quat& orient) { q_rocket = orient; } // set absolute orientation

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // functions used by sim                                                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    void update_dynamics(double current_time);
    void update_mass();
    void update_flight_controller(double current_time);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // functions used by flight controller                                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    void light_engine(); // should be used once per stage
    void cutoff_engine();
    void command_final_burn_fraction(double fraction); 
    bool advance_stage();
    void set_engine_orientation(Quat orientation);
    void rcs_on() { rcs_active = true; } // enable or disable rcs orientation correction
    void rcs_off() { rcs_active = false; }
    void rcs_apply_const_moment(Vec3 moment); // applies moment until changed
    void activate_detonation() { detonated = true; }

    // rocket owns the flight controller state
    Rocket(Rocket&&) = default;
    Rocket& operator=(Rocket&&) = default;
    Rocket(const Rocket&) = delete;
    Rocket& operator=(const Rocket&) = delete;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // private                                                                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // rocket static configuration                                                               //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // the dropped in flight controller
    INS ins;
    fc_sim_connector::State fc; // state fc_init
    fc_commands fc_cmd{}; // what it asked for on the current step

    std::unique_ptr<fc_vehicle> fc_veh;
    std::vector<fc_stage> fc_stages;

    bool fc_started = false;
    double fc_last_time = 0.0;

    // applies the commands sent by the FC through the API
    void apply_fc_commands();

    // topography
    renderer::EarthSurface* topo = &renderer::EarthSurface::Get();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // rocket static configuration                                                               //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    RocketProps props; // geometric properties of the rocket, set by the config
    std::string name;

    // flight data csv writer
    std::unique_ptr<DataExport> data_export;

    // initial launch geometry (origin, target, launch attitude, such things)
    RocketStartState start_state;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // dynamic state                                                                             //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int active_idx = 0;          // index of the currently active stage
    bool engine_locked = false;  // once cut off, the active stage's motor cannot be relit until staged away
    bool pending_cutoff = false; // a sub-step burn is finishing; thrust is zeroed at the start of the next step

    // mass properties
    double m_current = 0;        // current total mass (kg)
    double m_fuel_current = 0;   // current total fuel mass (kg)
    Vec3 I_body = {0, 0, 0};     // moments of inertia about the combined CoM, body frame
    double z_cm = 0;             // combined CoM along body +z, from the active stage's aft edge (m)
    double z_cp = 0;             // center of pressure along body +z, from the active stage's aft edge (m)
    Vec3 aero_torque = {0, 0, 0}; // aerodynamic moment about the combined CoM, body frame (N-m)

    // rcs system
    bool rcs_active = false;
    Vec3 applied_rcs_moment = {0, 0, 0}; // the moment that is applied to the body in addition to other forces due to RCS

    // kinematic state
    Vec3 r = {0, 0, 0};             // position (m)
    Vec3 v = {0, 0, 0};             // velocity (m/s)
    Vec3 a = {0, 0, 0};             // acceleration (m/s^2)
    Vec3 a_spec = {0, 0, 0};        // specific force in the body frame (m/s^2) (INS reads this)
    Vec3 w = {0, 0, 0};             // angular velocity (rad/s)
    Quat q_rocket = {1, 0, 0, 0};   // orientation of rocket nose relative to ECI (+z is nose)
    Quat q_engine = {1, 0, 0, 0};   // orientation of engine relative to rocket body
    double altitude = 0;

    // tracked meta properties for data analytics
    Vec3 drag_accel = {0, 0, 0};
    Vec3 grav_accel = {0, 0, 0};
    Vec3 thrust_accel = {0, 0, 0};
    double mach = 0;
    double dyn_pressure = 0;
    double aoa = 0;


    // accessors for the currently active stage
    Stage& active_stage() { return props.stages[active_idx]; }
    const Stage& active_stage() const { return props.stages[active_idx]; }
    int num_stages() const { return static_cast<int>(props.stages.size()); }

    // rocket explode button
    bool detonated = false;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // helper functions                                                                          //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // starting
    void set_start(double origin_latitude, double origin_longitude, double target_latitude, double target_longitude); // sets the starting and target position/attitude (only called from the constructor
    
    // kinematic helpers
    Vec3 engine_thrust_body(double thrust_scale) const;
    Vec3 net_body_torque(double thrust_scale) const; // engine + rcs torque about the combined CoM, body frame (constant across a step)
    void apply_ground_dynamics(const Vec3& I, double m_end, double dt);

    // applies translational acceleration components
    Vec3 translational_accel(double m_i, const Vec3& r_i, const Vec3& v_i, const Quat& q_i, const Vec3& w_i, const Vec3& thrust_body, const RocketProps& props); // gravity + drag + thrust, ECI
        Vec3 calc_drag_accel(const Vec3& r, const Vec3& v, const Quat& q, const Vec3& w, double mass, const RocketProps& props); // also sets aero_torque
        Vec3 calc_gravity_accel(const Vec3& r);


    // coordinate system conversion helpers
    Vec3 nose_direction_eci(const Quat& q) const;
    Vec3 lat_lon_to_ecef(double latitude_deg, double longitude_deg);

    // rocket state helpers
    double rocket_body_length() const; // nose to aft end of the remaining stack (m)
    bool is_rocket_on_ground(double com_dist_from_gnd); // snaps the rocket onto the surface if it is touching the ground
};

// standard atmosphere layers (air density/pressure, speed of sound, and dynamic viscosity mu at a given altitude above sea level)
void atmosphere(double altitude, double& air_density, double& air_pressure, double& speed_of_sound, double& mu);
