#pragma once

#include "render_types.hpp"
#include "rmath.hpp"
#include "../types.hpp"
#include <string>
#include <vector>

// Backend-neutral *domain* description of the scene's high-level objects. The
// renderer computes these each frame from the sim (it owns the precision-
// critical ECI->view shift).

namespace renderer {

// Static rocket dimensions (metres). Changes at staging; the backend rebuilds
// its meshes when these move.
struct RocketDims {
    double length;
    double cm_dist;
    double radius;
    double engine_dist;
    double nose_length;
};

// Everything needed to draw the rocket for one frame. Transforms are already in
// view space (the backend must not redo the ECI->view shift, or float jitter
// returns). The backend owns the actual geometry/materials.
struct RocketFrame {
    RocketDims dims;        // current stack dimensions
    RMat4 hull;            // body-local -> view
    RMat4 bell;            // bell-local (gimballed about the engine pivot) -> view
    bool  has_engine;      // false = stage has no engine
    bool  firing;          // plume visible
    float thrust;          // [0,1] plume intensity
    float flick;           // flame flicker multiplier
    float air;             // [0,1] atmospheric density (1 = sea level, 0 = vacuum)
    float heating;         // [0,1] aerodynamic heating (ablation glow); ascent/reentry
    RVec3 vel_dir;         // unit travel direction, view space (windward = normal . vel_dir)
    RVec3 nozzle;          // bell exit, view space (km)
    RVec3 exhaust_dir;     // unit direction the plume travels, view space

    // Detonation. When `detonated` is set the rocket has blown up: the backend
    // draws the explosion at `center` and no hull/plume, driving the animation
    // off `det_time`.
    bool  detonated = false;   // this rocket has detonated
    float det_time  = 0.0f;    // seconds elapsed since detonation began
    RVec3 center{};            // rocket position, view space (km) -- explosion anchor
};

// One rocket, as the HUD sees it.
struct HudRocket {
    std::string id;          // Greek ID, e.g. "α3"
    RVec3  view_pos;         // position, view space (km), for screen labels
    bool   detonated;
    float  thrust;           // [0,1] plume level (the renderer's firing estimate)
    double length;           // stack length (m); drops at staging
    double alt_km;           // above the sea-level sphere
    double vspeed;           // m/s along local up
};

// Everything the HUD shows, gathered once per frame after the 3D pass. A
// backend that draws its own HUD (RenderBackend::DrawHud) gets this in place
// of the renderer's default panels.
struct HudFrame {
    std::vector<HudRocket> rockets;   // every rocket, indexed like the sim's list
    int primary = -1;                 // selected rocket, -1 when none are live

    // Selected rocket, SI units unless noted.
    double met = 0, mass = 0, fuel = 0;
    Vec3   pos_km{};                  // ECEF
    double alt_km = 0, lat_deg = 0, lon_deg = 0;
    double speed = 0, vspeed = 0, accel = 0;
    double target_range_km = 0;       // great circle, ground point -> aim point
    double pitch_deg = 0;             // nose above the local horizon
    double gimbal_deg = 0;            // total nozzle deflection
    double gimbal_x_deg = 0, gimbal_y_deg = 0;   // deflection toward body +X / +Y
    Vec3   rates_dps{};               // body roll / pitch / yaw rates

    // Overlay toggles (number keys), in key order.
    struct Toggle { int key; const char* name; bool on; };
    std::vector<Toggle> toggles;
    bool show_telemetry = true;
    bool show_labels    = true;
};

// Everything needed to draw the Earth for one frame. The backend owns the mesh,
// texture(s), and shader; the renderer supplies the placement and lighting.
struct EarthFrame {
    RMat4 model;           // ECI-metre sphere -> view
    RVec3 sun_dir;         // unit view-space direction TO the sun
    RVec3 center;          // view-space sphere centre (km)
    RVec3 cam_pos;         // view-space camera position (km)
    // `center` again in double. At ~6400 km a float only resolves ~0.5 m, which
    // is visible between the ground and a rocket sitting on it; the terrain
    // anchors its chunks from this instead.
    Vec3  center_km;
};

}
