#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include "fc/inc/fc.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////
// entry points                                                                 //
///////////////////////////////////////////////////////////////////////////////////////////////

struct fc_state : FlightController {
    using FlightController::FlightController;
};

fc_state* fc_init(const fc_vehicle* vehicle) {
    return new fc_state(*vehicle);
}

void fc_update(fc_state* state, const fc_sensors* sensors, fc_commands* cmd) {
    state->flight_controller_process(*sensors, *cmd);
}

void fc_free(fc_state* state) {
    delete state;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// startup                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////////////

FlightController::FlightController(const fc_vehicle& vehicle) : veh(vehicle) {
    cs.stage = STANDBY;
}

FCInitState FlightController::create_target_trajectory(double lat_target, double long_target) {
    FCInitState out;

    // target's terrain height
    double target_radius = veh.r_target_ecef.mag();

    // convert to radians
    lat_target = lat_target * M_PI / 180.0;
    long_target = long_target * M_PI / 180.0;

    // derive lat and longitude from starting position of rocket on planet
    Vec3 p = veh.r_origin_eci;
    double radius = p.mag();
    double lat_origin  = asin(p.z / radius);
    double long_origin = atan2(p.y, p.x);

    // get ECI coordinates from lat and long
    out.r_origin = {
        radius * cos(lat_origin) * cos(long_origin),
        radius * cos(lat_origin) * sin(long_origin),
        radius * sin(lat_origin)
    };
    out.r_target_ecef = {
        target_radius * cos(lat_target) * cos(long_target),
        target_radius * cos(lat_target) * sin(long_target),
        target_radius * sin(lat_target)
    };

    Vec3 unit_vec_from_center = {
        .x = cos(lat_origin) * cos(long_origin),
        .y = cos(lat_origin) * sin(long_origin),
        .z = sin(lat_origin)
    };
    Vec3 up = {0, 0, 1};
    Vec3 rotation_axis = up.cross(unit_vec_from_center);
    double axis_norm = rotation_axis.mag();
    Vec3 rot_axis_u = axis_norm > 1e-12 ? rotation_axis / axis_norm : Vec3{1, 0, 0}; // any axis works at the poles
    double half_theta = acos(sin(lat_origin)) / 2;
    out.q_origin = {
        .w = cos(half_theta),
        .x = sin(half_theta) * rot_axis_u.x,
        .y = sin(half_theta) * rot_axis_u.y,
        .z = sin(half_theta) * rot_axis_u.z
    };

    // launch azimuth
    double delta_long = long_target - long_origin;
    double launch_azimuth = atan2(
        sin(delta_long) * cos(lat_target),
        cos(lat_origin) * sin(lat_target) - sin(lat_origin) * cos(lat_target) * cos(delta_long)
    );
    if (launch_azimuth < 0) launch_azimuth += 2.0 * M_PI;
    out.launch_asimuth = launch_azimuth;

    // burn time estimate for each stage
    out.stage_burn_time.reserve(num_stages());
    for (int i = 0; i < num_stages(); i++) {
        out.stage_burn_time.push_back(fc_stage_burn_time(&stage(i)));
    }

    // return all our calculated stuff yay!
    return out;
}

// init
void FlightController::init(double current_time) {
    cs = {};
    cs.stage = ARMED;
    cs.time = current_time;

    double tgt_lat = asin(veh.r_target_ecef.z / veh.r_target_ecef.mag()) * RAD_TO_DEG;
    double tgt_long = atan2(veh.r_target_ecef.y, veh.r_target_ecef.x) * RAD_TO_DEG;
    cs.is = create_target_trajectory(tgt_lat, tgt_long);

    cs.r = cs.is.r_origin; // set initial r to starting r
    cs.v = surface_velocity_eci(cs.r); // pad rotates with earth
    cs.att = cs.is.q_origin; // set initial orientation
    cs.target_att = cs.is.q_origin;
    countdown_start = current_time;

    cs.fuel.reserve(num_stages());
    for (int i = 0; i < num_stages(); i++) cs.fuel.push_back(stage(i).m_fuel);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// helper functions                                                                          //
///////////////////////////////////////////////////////////////////////////////////////////////

// aquires new data from the sim
void FlightController::pull_new_data(const fc_sensors& sensors) {
    cs.g = fc_gravity_j2(cs.r);
    cs.a_inertial = fc_q_rotate(cs.att, sensors.a_spec); // body -> ECI
    cs.a = cs.a_inertial + cs.g;
    cs.w = sensors.w;
    cs.dt = sensors.dt;
    cs.time = sensors.t;
}

// estimate state of the rocket in flight at the current moment
void FlightController::estimate_state() {
    // position and velocity 

    // acceleration is assumed to have been determined before this by reading from the INS

    // add delta v based on a to current v
    cs.v = cs.v + cs.a * cs.dt;

    // add delta r based on v to current r
    cs.r = cs.r + cs.v * cs.dt;

    // now for attitude estimation

    // integrate quaternion rate
    Quat q_dot = cs.att * Quat{0.0, cs.w.x, cs.w.y, cs.w.z};
    cs.att.w += 0.5 * q_dot.w * cs.dt;
    cs.att.x += 0.5 * q_dot.x * cs.dt;
    cs.att.y += 0.5 * q_dot.y * cs.dt;
    cs.att.z += 0.5 * q_dot.z * cs.dt;

    // normalize quaternion bc computer shit
    cs.att = cs.att.normalize();
}

// returns quaternion that points in direction of
Quat FlightController::quat_from_vec(Vec3 u) {
    u = u.unit();

    // quaternion of the shortest arc from nose to u
    Quat q = {
        .w = 1.0 + u.z,
        .x = -u.y,
        .y =  u.x,
        .z =  0.0,
    };
    return q.normalize();
}

// clamps a pointing direction to the angle of attack the airframe can hold at the current dynamic pressure
// without destabalizing
// @todo make this better
Vec3 FlightController::limit_aoa(Vec3 dir) const {
    Vec3 u = dir.unit();

    // airspeed relative to the atmosphere rotating with the earth
    Vec3 v_air = cs.v - surface_velocity_eci(cs.r);
    double speed = v_air.mag();
    if (speed < 1.0) return u;
    Vec3 v_hat = v_air * (1.0 / speed);

    // rough atmosphere model for dynamic pressure
    double h = cs.r.mag() - FC_EARTH_RADIUS;
    double rho = 1.225 * exp(-h / 8500.0);
    double q = 0.5 * rho * speed * speed;
    double aoa_max = MAX_Q_ALPHA / std::max(q, 1e-9);

    double aoa = acos(std::clamp(u.dot(v_hat), -1.0, 1.0));
    if (aoa <= aoa_max) return u;

    // point at the target direction, stopping at aoa_max
    Vec3 perp = u - v_hat * u.dot(v_hat);
    if (perp.mag() < 1e-9) return v_hat;
    return v_hat * cos(aoa_max) + perp.unit() * sin(aoa_max);
}

// sets the engine gimbal based on target orientation
Quat FlightController::set_new_engine_gimbal_quat() {
    if (cs.stage < STAGE_1 || cs.stage == FREE_FLIGHT) return {1, 0, 0, 0};

    Quat target = cs.target_att;
    Quat current = cs.att;

    // calculate error between the two quaternions
    Quat q_err = current.inverse() * target;
    if (q_err.w < 0) q_err = {-q_err.w, -q_err.x, -q_err.y, -q_err.z};

    // for small angles the vector part of the error quaternion is proportional 
    // to the rotational error about body axis. n is the roll, pitch, yaw error in radians
    Vec3 n = { .x = 2 * q_err.x, .y = 2 * q_err.y, .z = 2 * q_err.z };

    // calculate required torque using PD
    double wn = 5.0;
    double zeta = 0.7;
    Vec3 K_p = cs.I * (wn * wn);
    Vec3 K_d = cs.I * (2.0 * zeta * wn);
    Vec3 tau_req = n.vector_individual_multiply(K_p) - cs.w.vector_individual_multiply(K_d);

    // map torque to gimball command angles

    // active stage index for the current mission stage
    const fc_stage& s = stage(cs.active_stage);

    // CoM of the rocket
    double z_cm = cs.z_cm;

    // engine gimbal point
    double s_engine = s.tip_to_end_length - s.engine_distance;

    // moment arm from the engine gimbal point to the rocket CoM
    double moment_arm = z_cm - s_engine;

    // required pitch and yaw
    double pitch = -tau_req.x / (s.max_thrust * moment_arm);
    double yaw = -tau_req.y / (s.max_thrust * moment_arm);

    // determine nozzzle deflection from body in quaternion orientation from pitch and yaw commands
    Quat q_pitch = { cos(pitch / 2), sin(pitch / 2), 0.0, 0.0 };
    Quat q_yaw = { cos(yaw / 2), 0.0, sin(yaw / 2), 0.0 };
    return (q_pitch * q_yaw).normalize();
}

// figures out what moment needs to be applied by the RCS system to achieve the target orientation 
Vec3 FlightController::calculate_rcs_moments_to_achieve_target_orientation() {
    Quat target = cs.target_att;
    Quat current = cs.att;

    // calculate error between the two quaternions
    Quat q_err = current.inverse() * target;
    if (q_err.w < 0) q_err = {-q_err.w, -q_err.x, -q_err.y, -q_err.z};

    // for small angles the vector part of the error quaternion is proportional 
    // to the rotational error about body axis. n is the roll, pitch, yaw error in radians
    Vec3 n = { .x = 2 * q_err.x, .y = 2 * q_err.y, .z = 2 * q_err.z };

    // do PD for required orientation
    Vec3 I = cs.I;
    double wn = 0.15; // random placeholder
    double zeta = 0.7;
    Vec3 K_p = I * (wn * wn);
    Vec3 K_d = I * zeta * wn;
    Vec3 tau_req = n.vector_individual_multiply(K_p) - cs.w.vector_individual_multiply(K_d);
    //std::cout << "PD applying moment of " << tau_req.x << ", " << tau_req.y << ", " << tau_req.z << " n-m\n";
    return { tau_req.x, tau_req.y, tau_req.z };
}

// drains the propellant estimate for whatever burned over the last step
void FlightController::track_fuel() {
    if (cs.lit_stage < 0) return;

    double& f = cs.fuel[cs.lit_stage];
    f = std::max(0.0, f - fc_stage_max_mass_flow(&stage(cs.lit_stage)) * cs.throttle * cs.dt);

    // sub step burn only lasts one step
    if (cs.burn_ends) {
        cs.burn_ends = false;
        cs.lit_stage = -1;
        cs.throttle = 0.0;
    }
}

// fraction of a stage's propellant estimated to still be in the tank
double FlightController::fuel_fill(int i) const {
    return stage(i).m_fuel > 0 ? cs.fuel[i] / stage(i).m_fuel : 0.0;
}

// propellant column CoM from the stage tip as the tank drains
double FlightController::fuel_CoM(int i) const {
    return stage(i).fuel_CoM_dist + 0.5 * stage(i).fuel_length * (1.0 - fuel_fill(i));
}

// dry structure CoM from the stage tip
double FlightController::dry_CoM(int i) const {
    const fc_stage& st = stage(i);
    if (st.m_dry > 0) {
        return ((st.m_dry + st.m_fuel) * st.CoM_dist - st.m_fuel * st.fuel_CoM_dist) / st.m_dry;
    } else {
        return st.CoM_dist;
    }
}

// delta v a stage can give the rocket at a given point in time
double FlightController::stack_delta_v(int i) const {
    double m0 = veh.nosecone_mass;
    for (int j = i; j < num_stages(); j++) m0 += stage(j).m_dry + cs.fuel[j];
    double mf = m0 - cs.fuel[i];
    return mf > 0.0 ? fc_stage_exhaust_velocity(&stage(i)) * log(m0 / mf) : 0.0;
}

// integrates things to give a decent estimate of what the current moment of inertia of the rocket is
void FlightController::calculate_I() {
    Vec3 I = {0};
    double R2 = veh.radius * veh.radius;

    const int first_stage = cs.active_stage;
    const int stage_count = num_stages();

    // mass-weighted center of mass of the remaining stages
    double M_total = 0.0, m_CoM = 0.0, base = 0.0;
    for (int i = first_stage; i < stage_count; i++) {
        const fc_stage& st = stage(i);
        double L = st.tip_to_end_length;
        M_total += st.m_dry + cs.fuel[i];
        m_CoM += st.m_dry * (base + L - dry_CoM(i)) + cs.fuel[i] * (base + L - fuel_CoM(i));
        base += L;
    }

    // nosecone
    double m_nose = veh.nosecone_mass;
    double z_nose = base + veh.nosecone_length - veh.nosecone_com_distance;
    M_total += m_nose;
    m_CoM += m_nose * z_nose;

    double z_cm = m_CoM / M_total;
    cs.z_cm = z_cm;

    base = 0.0;
    for (int i = first_stage; i < stage_count; i++) {
        const fc_stage& s = stage(i); // selected stage
        double L = s.tip_to_end_length;
        double L_f = s.fuel_length * fuel_fill(i);

        // about z axis
        I.z += 0.5 * (s.m_dry + cs.fuel[i]) * R2;

        // about x and y axes parallel axis theorem from the structure and propellant centroids to z_cm
        double d_dry = (base + L - dry_CoM(i)) - z_cm;
        double d_fuel = (base + L - fuel_CoM(i)) - z_cm;
        double Ixy = s.m_dry * (L * L / 12.0 + R2 / 4) + s.m_dry * d_dry * d_dry
                   + cs.fuel[i] * (L_f * L_f / 12.0 + R2 / 4) + cs.fuel[i] * d_fuel * d_fuel;
        I.x += Ixy;
        I.y += Ixy;

        base += L;
    }

    // nosecone as a thin conical shell
    double h = veh.nosecone_length, d_nose = z_nose - z_cm;
    I.z += 0.5 * m_nose * R2;
    I.x += m_nose * (R2 / 4 + h * h / 18.0) + m_nose * d_nose * d_nose;
    I.y += m_nose * (R2 / 4 + h * h / 18.0) + m_nose * d_nose * d_nose;

    cs.I = I;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// flight controller loop                                                                    //
///////////////////////////////////////////////////////////////////////////////////////////////

// called once per sim time step
void FlightController::flight_controller_process(const fc_sensors& sensors, fc_commands& cmd) {

    ////////////////////////////
    // control inside         //
    ////////////////////////////

    // NOTE!! the staging functions inside the switch statement only ever set flags. the fc
    // commands all get issued in one block at the bottom so the behavior isnt hidden and we
    // dont accidentally do shit we dont want to

    if (cs.stage == STANDBY) init(sensors.t);

    // get latest data from the INS and advance the clock
    pull_new_data(sensors);

    // estimate state based on pulled data
    estimate_state();
    track_fuel();

    // manages setting a target attitude for the engine gimballing stuff depending on what the stage is
    switch (cs.stage) {
    case STANDBY:
        break;
    case ARMED:
        if (cs.time - countdown_start >= HOLD_DURATION) {
            cs.light_engine_flag = true;
            cs.stage_burn_time_start = cs.time;
            cs.stage = STAGE_1;
        }
        break;
    case STAGE_1: {
        s1_powered();
        break;
    }
    case STAGE_2:
        s2_powered();
        break;
    case PAYLOAD_DEPLOY:
        payload_deploy();
        break;
    case FREE_FLIGHT:
        free_flight();
        break;
    }

    // estimate current moment of inertia
    calculate_I();

    // check if the engine was supposed to be cut off
    if (cs.final_burn_flag) {
        cs.final_burn_flag = false;
        // fractional burn (burn that is sub step to time step)
        cmd.cutoff = 1;
        cmd.cutoff_fraction = cs.final_burn_fraction;
        cs.throttle = std::clamp(cs.final_burn_fraction, 0.0, 1.0);
        cs.burn_ends = true;
    }
    else if (cs.cutoff_engine_flag) {
        cs.cutoff_engine_flag = false;
        cmd.cutoff = 1;
        cs.lit_stage = -1;
        cs.throttle = 0.0;
    }

    // check if the stage was supposed to be separated (clears the engine lock for the fresh stage)
    if (cs.separate_stage_flag) {
        cs.separate_stage_flag = false;
        cmd.separate = 1;
        if (cs.active_stage + 1 < num_stages()) cs.active_stage++;
        cs.lit_stage = -1;
        cs.throttle = 0.0;
    }

    // check if engine was supposed to be lit
    if (cs.light_engine_flag) {
        cs.light_engine_flag = false;
        cmd.light = 1;
        if (cs.lit_stage != cs.active_stage) {
            cs.lit_stage = cs.active_stage;
            cs.throttle = 1.0;
        }
    }

    // check if rocket was supposed to be detonated
    if (cs.detonate_flag) {
        cmd.detonate = 1;
    }

    // send targeting commands to the engine gimbal system based on target attitude in cs
    cmd.gimbal = set_new_engine_gimbal_quat();

    // send orientation change commands to the rcs thruster system if it is active
    if (cs.rcs_activated_flag) {
        cmd.rcs_on = 1;
        cmd.rcs_moment = calculate_rcs_moments_to_achieve_target_orientation();
    }
    else {
        cmd.rcs_on = 0;
    }
}
