#include "fc/inc/fc.hpp"
#include <iostream>
#include <algorithm>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// S1 GUIDANCE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// rough init tof guess
constexpr double NOMINAL_TOF = 1500.0; // seconds

void FlightController::s1_powered() {
    double dt = cs.time - cs.stage_burn_time_start;

    // aim at where the target will roughly be at arrival
    Vec3 r_target = target_eci_at_time_of_arrival(cs.time + NOMINAL_TOF);

    ///////////////////////////////////////////////////////////////////////////
    // do initial turn towards target before gravity turn starts
    ///////////////////////////////////////////////////////////////////////////
    constexpr double init_turn_len = 20.0; // seconds
    constexpr double init_tilt_angle = 18.0 * DEG_TO_RAD;

    if (dt < init_turn_len) {
        Vec3 up = cs.r.unit();

        // direction towards the target in the plane perpendicular to the rockets starting up point
        Vec3 downrange_from_origin = r_target - cs.is.r_origin;
        downrange_from_origin = (downrange_from_origin - up * downrange_from_origin.dot(up)).unit();

        // turned slightly
        Vec3 tilt = up * cos(init_tilt_angle) + downrange_from_origin * sin(init_tilt_angle);
        cs.target_att = quat_from_vec(limit_aoa(tilt));
    }

    ///////////////////////////////////////////////////////////////////////////
    // after initial turn, start gravity turn (follow v vec)
    ///////////////////////////////////////////////////////////////////////////
    else {
        // align with plane connecting to target and blend with current velocity vector relative to the air
        Vec3 n = cs.is.r_origin.cross(r_target).unit();
        Vec3 v_air = cs.v - surface_velocity_eci(cs.r);
        Vec3 v_in_target_plane = v_air - n * v_air.dot(n);
        cs.target_att = quat_from_vec(limit_aoa(v_in_target_plane));
    }

    ///////////////////////////////////////////////////////////////////////////
    // separate once the stage has burned out
    ///////////////////////////////////////////////////////////////////////////
    // note:: you should add stage separation additionally take into account
    // vibration and force checking to see if the engines are on thru the IMU
    // blended with this estimate!
    if (cs.time - cs.stage_burn_time_start > cs.is.stage_burn_time[0] + 1) {
        cs.separate_stage_flag = true;
        cs.light_engine_flag = true;
        cs.stage = STAGE_2;
        cs.stage_burn_time_start = cs.time;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// S2 GUIDANCE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// stumpff functions helper
static double Cz(double z) {
    if (z >  1e-6) return (1.0 - cos(sqrt(z))) / z;
    if (z < -1e-6) return (cosh(sqrt(-z)) - 1.0) / (-z);
    return 0.5;
}
static double Sz(double z) {
    if (z >  1e-6) { 
        double s = sqrt(z);
        return (s - sin(s)) / (s*s*s);
    }
    if (z < -1e-6) { 
        double s = sqrt(-z);
        return (sinh(s) - s) / (s*s*s);
    }
    return 1.0 / 6.0;
}

// helper for y(z) function in lambert
static double lambert_y(double z, double r_norm, double r_t_norm, double A) {
    return r_norm + r_t_norm + A * (z * Sz(z) - 1.0) / sqrt(Cz(z));
}

// time of flight error at a given z
static double lambert_F(double z, double r_norm, double r_t_norm, double A, double tff) {
    double y = lambert_y(z, r_norm, r_t_norm, A);
    return pow(y / Cz(z), 1.5) * Sz(z) + A * sqrt(y) - sqrt(GM_EARTH) * tff;
}

// lambert problem (determine velocity vector that reaches target in time tff given current position)
static Vec3 lambert(Vec3 r, Vec3 r_t, double tff) {
    double r_norm = r.mag();
    double r_t_norm = r_t.mag();

    // transfer angle
    double d_theta = acos(r.dot(r_t) / (r_norm * r_t_norm));

    // some constant
    double A = sin(d_theta) * sqrt( (r_norm * r_t_norm) / (1 - cos(d_theta)) );

    // iterate z until the trajectory's time of flight is the same as the input time of flight
    double z = 0.0; // initial number guess
    while (lambert_y(z, r_norm, r_t_norm, A) < 0.0) z += 0.1; // bump up z until y(z) is good enough

    // now iterate on z to find the best one using secant search
    double z0 = 0.0;
    double z1 = 0.1; // start with initial guesses that are kinda close to one another

    double F0 = lambert_F(z0, r_norm, r_t_norm, A, tff);
    for (int i = 0; i < 60; i++) {
        double F1 = lambert_F(z1, r_norm, r_t_norm, A, tff);

        // step
        double z_next = z1 - F1 * (z1 - z0) / (F1 - F0);

        // check if it converges
        if (fabs(z_next - z1) < 1e-8) {
            z = z_next;
            break;
        }

        // shift values for next iteration
        z0 = z1;
        F0 = F1;
        z1 = z_next;
    }

    // now we determine the lagrange multipliers as part of the analytical version of the lambert problem to find our required velocity vector
    double y = lambert_y(z, r_norm, r_t_norm, A);
    double f = 1 - y / r_norm;
    double g = A * sqrt(y / GM_EARTH);

    return (r_t - r * f) * (1 / g); // optimal taret velocity
}

// velocity needed to hit where the target will be after tof of earth spin
Vec3 FlightController::v_req_for_tof(double tof) const {
    return lambert(cs.r, target_eci_at_time_of_arrival(cs.time + tof), tof);
}

// delta v needed to get onto the tof trajectory
double FlightController::dv_for_tof(double tof) const {
    return (v_req_for_tof(tof) - cs.v).mag();
}

void FlightController::s2_powered() {
    cs.s2_lambert_counter++;
    double burn_time = cs.time - cs.stage_burn_time_start;

    Vec3& v_req = cs.s2_v_req; // optimal required velocity

    ////////////////////////////////////////////////////////////
    // determine the minimum rquired velocity
    ////////////////////////////////////////////////////////////
    if (cs.s2_lambert_counter >= cs.s2_lambert_counter_reset_num) {
        cs.s2_lambert_counter = 0; // reset counter

        ////////////////////////////////////////////////////////////////////////////////////////
        // use golden section search to find the optimal required velocity vector using lambert problem
        ////////////////////////////////////////////////////////////////////////////////////////
        double rho = 1 / ( (1 + sqrt(5) ) / 2);

        double a = 300; // lower bound of guess
        double b = 2700; // no reasonable flight should take longer than this right?

        // pick two points based on upper and lower bounds
        double x1 = b - rho * (b - a);
        double x2 = a + rho * (b - a);

        // do the lambert problem on these two points (smallest change in velocity)
        double f1 = dv_for_tof(x1);
        double f2 = dv_for_tof(x2);

        // iterate on a and b to find the minimum delta v point
        while ((b - a) > 0.1) {
            if (f1 < f2) {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - rho * (b - a);
                f1 = dv_for_tof(x1);
            }
            else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + rho * (b - a);
                f2 = dv_for_tof(x2);
            }
        }

        v_req = v_req_for_tof(0.5 * (a + b));
    }

    // velocity to be gained as a target
    Vec3 v_gain = v_req - cs.v;

    // stop the engine if the V is within proper cutoff range OR if the fuel is out
    // the cutoff should be tied to the max possible delta V of the 3rd stage engine
    double next_delta_v = stack_delta_v(2); // the 3rd stage also has to push the nosecone
    double estimated_s2_burn_time = fc_stage_burn_time(&stage(1));
    if (v_gain.mag() < next_delta_v - 0.5 * next_delta_v || burn_time > estimated_s2_burn_time) { // within 20% for margin of error, or motor depleted
        // stop engine to stop overshoot
        cs.cutoff_engine_flag = true;
        std::cout << "ENGINE_CUTOFF at " << v_gain.mag() << "m/s " << "with next stage having delta v: " << next_delta_v << "m/s\n";
        
        // switch stage
        cs.stage = PAYLOAD_DEPLOY;

        // set stage separation
        cs.separate_stage_flag = true;

        return;
    }

    // set attitude to new target
    cs.target_att = quat_from_vec(limit_aoa(v_gain));

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// S3 GUIDANCE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// simply operates a modified version of stage 2 that improves trajectory accuracy with small engine
void FlightController::payload_deploy() {
    if (cs.payload_deploy_cutoff_done) return;

    cs.s3_lambert_counter++;

    // initially, set rcs flag to true so we can orient rocket in right direction
    cs.rcs_activated_flag = true;

    ////////////////////////////////////////////////////////////////////////////////////////
    // calculate error in direction to know when its safe to light engine (otherwise use RCS)
    ////////////////////////////////////////////////////////////////////////////////////////
    Quat target = cs.target_att;
    Quat current = cs.att;

    // calculate error between the two quaternions
    Quat q_err = current.inverse() * target;
    if (q_err.w < 0) q_err = {-q_err.w, -q_err.x, -q_err.y, -q_err.z};

    // for small angles the vector part of the error quaternion is proportional
    // to the rotational error about body axis. n is the roll, pitch, yaw error in radians
    Vec3 n = { .x = 2 * q_err.x, .y = 2 * q_err.y, .z = 2 * q_err.z };

    // burn is ready once the attitude error has gotten low enough so engins is not lit while the rocket is still wobbling around
    constexpr double attitude_tolerance = 0.5 * DEG_TO_RAD; // rad
    bool burn_ready = n.mag() < attitude_tolerance;
    if (burn_ready) { cs.light_engine_flag = true; }

    ////////////////////////////////////////////////////////////////////////////////////////
    // use golden section search to find the optimal required velocity vector using lambert problem 
    ////////////////////////////////////////////////////////////////////////////////////////
    Vec3& v_req = cs.s3_v_req; // optimal required velocity

    // only run the search every few steps
    if (cs.s3_lambert_counter >= cs.s3_lambert_counter_reset_num) {
        cs.s3_lambert_counter = 0; // reset counter

        double rho = 1 / ( (1 + sqrt(5) ) / 2);

        double a = 300; // lower bound of guess
        double b = 2700; // no reasonable flight should take longer than this right?

        // pick two points based on upper and lower bounds
        double x1 = b - rho * (b - a);
        double x2 = a + rho * (b - a);

        // do the lambert problem on these two points
        double f1 = dv_for_tof(x1);
        double f2 = dv_for_tof(x2);

        // iterate on a and b to find the minimum delta v point
        while ((b - a) > 0.1) {
            if (f1 < f2) {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - rho * (b - a);
                f1 = dv_for_tof(x1);
            }
            else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + rho * (b - a);
                f2 = dv_for_tof(x2);
            }
        }

        double tof = 0.5 * (a + b);
        v_req = v_req_for_tof(tof);
        cs.s3_t_arrival = cs.time + tof;
    }

    // velocity to be gained as a target
    Vec3 v_gain = v_req - cs.v;

    // cut off when the remaining v to gain is smaller than the delta v the engine will add in the next step
    double dv_next_step = cs.a_inertial.mag() * cs.dt;
    bool burning = cs.lit_stage == cs.active_stage && dv_next_step > 0.01;

    // remaining delta v projected onto the thrust axis
    double v_needed = v_gain.dot(cs.a_inertial.unit());

    // once a single full step would meet or pass the target, burn only the needed
    // fraction of this step and then cut off
    // we do this because there is a coarse time step and the burn cuttoff always falls somewhere between
    // time step increments, so this allows us to burn in between that by interpolation (epic sigma0)
    if (burning && v_needed < dv_next_step) {
        double frac = std::clamp(v_needed / dv_next_step, 0.0, 1.0);
        cs.final_burn_fraction = frac;
        cs.final_burn_flag = true;
        cs.rcs_activated_flag = false;
        cs.payload_deploy_cutoff_done = true;
        cs.stage = FREE_FLIGHT;
        std::cout << "ENGINE_CUTOFF (frac " << frac << ") residual " << (v_needed - frac * dv_next_step) << " m/s\n";
        return;
    }

    // make rocket not demolish its trajectory when its close to v cutoff
    if (!burning || v_gain.mag() > 3.0 * dv_next_step) {
        cs.target_att = quat_from_vec(v_gain.unit());
    }
}

void FlightController::free_flight() {
    // determine altitude relative to target
    Vec3 target_eci = target_eci_at_time_of_arrival(cs.time);
    double dist = (cs.r - target_eci).mag();

    if (dist < 1000.0) {
        cs.detonate_flag = true;
    }
}