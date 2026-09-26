#include "models.hpp"
#include "../geometry.hpp"
#include "../../constants.hpp"
#include <vector>
#include <cmath>

namespace renderer {

namespace {
constexpr int kSides = 32;   // hull/bell facet count

// An open, fluted (ridged) engine bell: a single-walled surface of revolution
// from throat (zT) to exit (zE), flaring outward, with longitudinal ridges.
// Drawn unculled so the concave interior is visible (the actual opening).
void appendBell(Mesh& m, float zT, float zE, float rT, float rE,
                int sides, int ribs, float amp, RColor col) {
    const int M = 14;
    uint32_t base = (uint32_t)m.verts.size();
    for (int i = 0; i <= M; ++i) {
        float t  = (float)i / M;
        float z  = zT + (zE - zT) * t;
        float Ro = rT + (rE - rT) * powf(t, 1.5f);
        float dRodz = (zE != zT) ? ((rE - rT) * 1.5f * powf(fmaxf(t,1e-4f), 0.5f)) / (zE - zT) : 0.0f;
        for (int j = 0; j <= sides; ++j) {
            float th = (float)(TAU * j / sides), cs = cosf(th), sn = sinf(th);
            float mod = 1.0f + amp * cosf(ribs * th);
            float R   = Ro * mod;
            float dRdth = Ro * (-amp * ribs * sinf(ribs * th));
            RVec3 dPt = { dRdth*cs - R*sn, dRdth*sn + R*cs, 0 };
            RVec3 dPz = { dRodz*mod*cs, dRodz*mod*sn, 1 };
            RVec3 n = rmath::normalize(rmath::cross(dPz, dPt));
            if (n.x*cs + n.y*sn < 0.0f) n = { -n.x, -n.y, -n.z };   // outward
            m.verts.push_back({ {R*cs, R*sn, z}, n, 0, 0, col });
        }
    }
    const int stride = sides + 1;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < sides; ++j) {
            uint32_t a = base + i*stride + j, b = a + 1;
            uint32_t c = base + (i+1)*stride + j, d = c + 1;
            m.idx.push_back(a); m.idx.push_back(b); m.idx.push_back(d);
            m.idx.push_back(a); m.idx.push_back(d); m.idx.push_back(c);
        }
}

bool sameDims(const RocketDims& a, const RocketDims& b) {
    return a.length == b.length &&
           a.radius == b.radius && a.engine_dist == b.engine_dist &&
           a.nose_length == b.nose_length;
}

} // namespace

void RocketModel::Ensure(RenderBackend& b, const RocketDims& dims) {
    if (cone_ == 0) {   // first use: build the dimension-independent meshes once
        cone_   = b.CreateMesh(geom::buildCone(kSides));
        sphere_ = b.CreateMesh(geom::buildSphere(1.0f, 16, 24));
    }
    // Already the selected hull/bell? (fast path for a homogeneous fleet)
    if (hull_ != 0 && sameDims(dims_, dims)) return;

    // Reuse a cached build for these dims if we've seen them, else build once.
    for (const HullBell& e : cache_) {
        if (sameDims(e.dims, dims)) { hull_ = e.hull; bell_ = e.bell; dims_ = dims; return; }
    }
    HullBell e = buildHullBell(b, dims);
    cache_.push_back(e);
    hull_ = e.hull; bell_ = e.bell; dims_ = dims;
}

RocketModel::HullBell RocketModel::buildHullBell(RenderBackend& b, const RocketDims& d) const {
    const float L      = (float)d.length;
    const float radius = (float)d.radius;

    const RColor kBody   = { 226, 229, 235, 255 };  // brushed silver
    const RColor kNose   = { 196,  58,  58, 255 };  // red cap
    const RColor kCollar = {  74,  80,  92, 255 };  // gunmetal collar
    const RColor kBell   = {  46,  47,  53, 255 };  // dark nozzle

    // Body frame: +Z = nose, nose tip fixed at the origin (cm-independent so the
    // mesh caches per stage). Draw slides it by cm_dist + nose_length to seat the
    // CoM at st.r. The nose is the sim's nosecone, on top of the stack (L doesn't
    // include it), so the tip stays the same size across staging.
    const float nose_z   = 0.0f;
    const float noseLen  = (float)d.nose_length;
    const float nose_base = nose_z - noseLen;
    const float tail_z   = nose_base - L;

    // PBR params (metallic, roughness) are packed into vertex UVs per part.
    auto setMR = [](Mesh& m, size_t from, float metal, float rough) {
        for (size_t i = from; i < m.verts.size(); ++i) { m.verts[i].u = metal; m.verts[i].v = rough; }
    };

    Mesh hull;
    size_t k;

    // Main body — glossy white paint: mostly diffuse with only a mild reflection
    // (low metallic so the white albedo dominates instead of acting as a mirror).
    k = hull.verts.size();
    geom::appendFrustum(hull, tail_z, nose_base, radius, radius, kSides, kBody);
    setMR(hull, k, 0.12f, 0.30f);

    // Thin collar where the nose meets the body.
    k = hull.verts.size();
    geom::appendFrustum(hull, nose_base - radius * 0.18f, nose_base,
                        radius * 1.03f, radius, kSides, kCollar, false, false);
    setMR(hull, k, 0.80f, 0.45f);

    // Smoothed (elliptical) nose — painted dielectric.
    std::vector<RVec3> prof;
    const int NS = 22;
    for (int i = 0; i <= NS; ++i) {
        float t = (float)i / NS;                       // 0 base .. 1 tip
        float z = nose_base + t * noseLen;
        float r = radius * sqrtf(fmaxf(0.0f, 1.0f - t * t));
        prof.push_back({ z, r, 0 });
    }
    k = hull.verts.size();
    geom::appendRevolve(hull, prof, kSides, kNose, /*capBase=*/true);
    setMR(hull, k, 0.05f, 0.5f);

    MeshHandle hullH = b.CreateMesh(hull);

    // Engine bell: open, fluted nozzle in the gimbal-pivot frame (throat at z=0,
    // exit 1.6 m down). Throat narrow, flaring to a wide ridged exit. Dark, rough metal.
    Mesh bell;
    appendBell(bell, 0.0f, -1.6f, radius * 0.45f, radius * 1.25f, kSides, 20, 0.05f, kBell);
    setMR(bell, 0, 0.55f, 0.55f);
    MeshHandle bellH = b.CreateMesh(bell);

    return { d, hullH, bellH };
}

void RocketModel::Draw(RenderBackend& b, const RocketFrame& f) const {
    // Destroyed: the explosion replaces the intact hull + plume entirely.
    if (f.detonated) { drawDetonation(b, f); return; }

    // The hull mesh is built nose-at-origin; cm_dist is measured from the top of
    // the stack, below the nosecone, so slide it down by both to land the CoM at
    // st.r (which f.hull maps to the mesh origin).
    RMat4 hullM = rmath::mul(f.hull, rmath::translate({0, 0, (float)(f.dims.cm_dist + f.dims.nose_length)}));

    Material body;  body.lit = true; body.cull = false;   // unculled: open collar + nozzle shapes
    b.DrawModel(hull_, hullM, body);
    Material bell;  bell.lit = true; bell.cull = false;    // unculled: concave nozzle interior
    if (f.has_engine) b.DrawModel(bell_, f.bell, bell);   // engineless stage (e.g. nosecone): no bell

    // Ablation sheath: re-draw the hull slightly inflated as an additive glow
    // shell, so the plasma follows the rocket's shape (windward incandescence +
    // a bulge ahead of the nose) instead of a sphere. Intensity comes from
    // u_heat in the shader; shown on ascent and (brightly) on reentry.
    if (f.heating > 0.03f) {
        Material heat; heat.emissive = true; heat.depth_write = false;
        b.DrawModel(hull_, rmath::mul(hullM, rmath::scale(1.05f)), heat);
    }

    if (!f.firing) return;

    // Exhaust plume: nested additive cones streaming out the nozzle along the
    // exhaust direction, depth-write off so the layers glow through each other.
    const double radius = dims_.radius;
    const double Lm     = dims_.length * 0.8 * f.thrust * f.flick;  // plume length, metres
    auto cone = [&](double r0_m, double len_m, RColor c) {
        c.a = (unsigned char)(c.a * f.thrust);
        RVec3 apex = { f.nozzle.x + f.exhaust_dir.x * (float)(len_m * M_TO_KM),
                       f.nozzle.y + f.exhaust_dir.y * (float)(len_m * M_TO_KM),
                       f.nozzle.z + f.exhaust_dir.z * (float)(len_m * M_TO_KM) };
        Material m; m.color = c; m.blend = BlendMode::Additive; m.depth_write = false;
        b.DrawModel(cone_, rmath::orientCone(f.nozzle, apex, (float)(r0_m * M_TO_KM)), m);
    };
    cone(radius * 1.5,  Lm * 1.25, { 180,  70, 20,  90 });   // outer haze
    cone(radius * 1.0,  Lm,        { 255, 140, 40, 140 });   // flame body
    cone(radius * 0.55, Lm * 0.6,  { 255, 235, 180, 210 });  // white-hot core
    cone(radius * 0.30, Lm * 0.32, { 205, 225, 255, 220 });  // blue-white inner jet

    // Mach diamonds: periodic shock nodes along the supersonic exhaust. They form
    // only in atmosphere (over/under-expanded nozzle), so fade with altitude
    // (f.air) and down the plume. In vacuum (f.air -> 0) they vanish.
    if (f.air > 0.12f) {
        const int   nD      = 5;
        const float spacing = (float)(radius * 1.7) * (0.15f + 0.5f * (1.0f - f.air));  // metres
        for (int k = 1; k <= nD; ++k) {
            float dist = spacing * k;
            if (dist > Lm * 0.95) break;
            float fade = 1.0f - (float)(k - 1) / nD;
            float br   = f.air * f.thrust * fade;
            RVec3 c = { f.nozzle.x + f.exhaust_dir.x * (float)(dist * M_TO_KM),
                        f.nozzle.y + f.exhaust_dir.y * (float)(dist * M_TO_KM),
                        f.nozzle.z + f.exhaust_dir.z * (float)(dist * M_TO_KM) };
            float sz = (float)(radius * 0.45 * M_TO_KM) * (0.6f + 0.4f * fade);
            Material d;
            d.color = { 215, 230, 255, (unsigned char)(220.0f * br) };
            d.blend = BlendMode::Additive; d.depth_write = false;
            b.DrawModel(sphere_, rmath::placeSphere(c, sz), d);
        }
    }

    // Nozzle glow.
    const float radiusW = (float)(radius * M_TO_KM);
    Material g;
    g.color = { 255, 190, 90, (unsigned char)(150 * f.thrust) };
    g.blend = BlendMode::Additive; g.depth_write = false;
    b.DrawModel(sphere_, rmath::placeSphere(f.nozzle, radiusW * (0.7f + 0.2f * f.flick)), g);
}

void RocketModel::drawDetonation(RenderBackend& b, const RocketFrame& f) const {
    const float a      = f.det_time;    // seconds since detonation
    const float Dfire  = 1.5f;          // fireball lifetime
    const float Dtotal = 4.0f;          // incl. lingering smoke
    if (a > Dtotal) return;             // fully dissipated: the rocket is gone

    auto clamp01 = [](float x) { return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x); };

    // Blast scale: a few tens of metres, floored so small rockets still read. All
    // layers are additive glow spheres about the rocket's view-space position,
    // depth-tested but not depth-writing (so the Earth occludes them and the layers
    // blend through each other), exactly like the plume.
    const float Rmax  = (float)(fmax(dims_.radius * 25.0, 30.0) * M_TO_KM);  // km
    const RVec3 c     = f.center;
    const float flick = f.flick;        // reuse the plume's flame flicker for turbulence

    auto glow = [&](float radiusKm, RColor rgb, float alpha) {
        alpha = clamp01(alpha);
        if (alpha <= 0.002f || radiusKm <= 0.0f) return;
        Material m;
        m.color = { rgb.r, rgb.g, rgb.b, (unsigned char)(alpha * 255.0f) };
        m.blend = BlendMode::Additive;
        m.depth_write = false;
        b.DrawModel(sphere_, rmath::placeSphere(c, radiusKm), m);
    };

    // Flash: brief blinding light at t=0.
    if (a < 0.12f) {
        float p = a / 0.12f;
        glow(Rmax * (0.5f + 0.6f * p), { 255, 250, 235, 255 }, (1.0f - p) * 0.9f);
    }

    // Fireball: nested layers, ease-out expansion, fading over Dfire.
    if (a < Dfire) {
        float p    = a / Dfire;
        float R    = Rmax * (1.0f - (1.0f - p) * (1.0f - p));   // fast then slowing
        float env  = 1.0f - p;                                  // overall fade envelope
        float turb = 0.9f + 0.15f * flick;
        glow(R * 1.15f * turb, { 120,  40,  10, 255 }, env * 0.35f);            // smoky outer
        glow(R * 0.85f * turb, { 255, 120,  30, 255 }, env * 0.60f);            // flame body
        glow(R * 0.50f,        { 255, 230, 180, 255 }, (1.0f - p * p) * 0.85f); // hot core
        if (p < 0.3f)                                                           // blue-white centre
            glow(R * 0.28f,    { 210, 225, 255, 255 }, (1.0f - p / 0.3f) * 0.7f);
    }

    // Shockwave: thin fast front, gone by 0.6 s.
    if (a < 0.6f) {
        float p = a / 0.6f;
        glow(Rmax * (0.4f + 1.8f * p), { 200, 220, 255, 255 }, (1.0f - p) * 0.25f);
    }

    // Lingering smoke: a dim warm glow that slowly grows and fades out by Dtotal.
    const float s0 = Dfire * 0.6f;      // starts overlapping the fireball tail
    if (a > s0) {
        float q = clamp01((a - s0) / (Dtotal - s0));
        glow(Rmax * (1.0f + 0.6f * q), { 60, 35, 20, 255 }, (1.0f - q) * 0.30f);
    }
}

}
