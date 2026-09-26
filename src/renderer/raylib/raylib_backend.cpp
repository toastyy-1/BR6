#include "raylib_backend.hpp"
#include "theme.hpp"
#include <raymath.h>
#include <rlgl.h>
#include <vector>

namespace renderer {

namespace {

inline Vector3 toRl(const RVec3& v) { return { v.x, v.y, v.z }; }
inline Color   toRl(RColor c)       { return { c.r, c.g, c.b, c.a }; }

// raylib's Matrix field m{k} is column-major element k, like RMat4's m[k]; the
// structs' memory orders differ, so copy by name (not memcpy).
RMat4 fromRl(const Matrix& r) {
    return RMat4{{ r.m0,  r.m1,  r.m2,  r.m3,  r.m4,  r.m5,  r.m6,  r.m7,
                   r.m8,  r.m9,  r.m10, r.m11, r.m12, r.m13, r.m14, r.m15 }};
}

} // namespace

void RaylibBackend::Init(int width, int height, const char* title) {
    // MSAA keeps the wireframe lines smooth. 60 fps is plenty to watch a flight
    // and keeps the renderer's share of the machine small on high-refresh
    // displays; vsync stops tearing.
    // raylib logs an INFO line for every buffer it frees. With a full terrain
    // cache (near the ground) that's ~10k lines at shutdown, and the window hangs
    // unresponsive while a slow terminal drains them. Warnings and errors only.
    SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT | FLAG_VSYNC_HINT);
    InitWindow(width, height, title);
    SetWindowMonitor(0);
    SetTargetFPS(60);

    wire_.Init();
    earth_.Init();
    coast_.Init();
    hud_.Init();
}

void RaylibBackend::Shutdown() {
    earth_.Shutdown();
    coast_.Shutdown();
    for (wire::GpuMesh& m : meshes_) wire::Destroy(m);
    for (::Texture2D& t : textures_) UnloadTexture(t);
    wire_.Shutdown();
    hud_.Shutdown();
    CloseWindow();
}

bool RaylibBackend::ShouldClose() const { return WindowShouldClose(); }

FrameInput RaylibBackend::PollInput() {
    Vector2 d = GetMouseDelta();
    float mx = (IsKeyDown(KEY_D) ? 1.0f : 0.0f) - (IsKeyDown(KEY_A) ? 1.0f : 0.0f);
    float my = (IsKeyDown(KEY_E) ? 1.0f : 0.0f) - (IsKeyDown(KEY_Q) ? 1.0f : 0.0f);
    float mz = (IsKeyDown(KEY_W) ? 1.0f : 0.0f) - (IsKeyDown(KEY_S) ? 1.0f : 0.0f);
    bool boost = IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT);
    bool recenter = IsKeyDown(KEY_F);
    // Tab cycles the primary rocket (edge-triggered); Shift+Tab goes backwards.
    bool tab  = IsKeyPressed(KEY_TAB);
    bool next = tab && !boost;
    bool prev = tab &&  boost;
    // Number keys 1-9 toggle overlays (edge-triggered). KEY_ONE..KEY_NINE are
    // consecutive, so map the first that was pressed this frame.
    int toggle = 0;
    for (int k = KEY_ONE; k <= KEY_NINE; ++k) if (IsKeyPressed(k)) { toggle = k - KEY_ONE + 1; break; }
    return FrameInput {
        d.x, d.y,
        GetMouseWheelMove(),
        IsMouseButtonDown(MOUSE_BUTTON_LEFT),
        mx, my, mz, boost, recenter,
        next, prev, toggle,
    };
}

float  RaylibBackend::FrameTime() const { return GetFrameTime(); }
double RaylibBackend::Time() const      { return GetTime(); }

TextureHandle RaylibBackend::LoadTexture(const char* path) {
    ::Texture2D t = ::LoadTexture(path);
    GenTextureMipmaps(&t);
    SetTextureFilter(t, TEXTURE_FILTER_TRILINEAR);
    textures_.push_back(t);
    return (TextureHandle)textures_.size();
}

MeshHandle RaylibBackend::CreateMesh(const Mesh& m) {
    meshes_.push_back(wire::Upload(m.verts, m.idx, wire::FeatureEdges(m)));
    return (MeshHandle)meshes_.size();
}

void RaylibBackend::DestroyMesh(MeshHandle h) {
    if (h == 0 || h > meshes_.size()) return;
    wire::Destroy(meshes_[h - 1]);   // tombstone (draws skip it); handle is not reused
}

void RaylibBackend::BeginFrame(RColor clear) {
    clear_ = clear;
    labels_.clear();
    BeginDrawing();
    ClearBackground(toRl(clear));
}

void RaylibBackend::EndFrame() {
    EndDrawing();
    wire_.NextFrame();
}

void RaylibBackend::SetClipPlanes(float near_plane, float far_plane) {
    rlSetClipPlanes(near_plane, far_plane);
}

void RaylibBackend::Begin3D(const RCamera& cam) {
    cam_   = cam;
    cam3d_ = Camera3D {
        toRl(cam.position),
        toRl(cam.target),
        toRl(cam.up),
        cam.fovy,
        CAMERA_PERSPECTIVE,
    };
    BeginMode3D(cam3d_);
    // BeginMode3D has just loaded this camera's projection and view into rlgl.
    wire_.Begin(rmath::mul(fromRl(rlGetMatrixProjection()), fromRl(rlGetMatrixModelview())), clear_);
}

void RaylibBackend::End3D() {
    wire_.End();
    EndMode3D();
}

ScreenPoint RaylibBackend::WorldToScreen(const RVec3& viewPos) const {
    // Cull points behind the camera: GetWorldToScreen still returns coordinates
    // for them, but mirrored/nonsensical, so gate on the forward half-space.
    Vector3 fwd = Vector3Normalize(Vector3Subtract(cam3d_.target, cam3d_.position));
    Vector3 rel = Vector3Subtract(toRl(viewPos), cam3d_.position);
    if (Vector3DotProduct(rel, fwd) <= 0.0f) return { 0, 0, false };
    Vector2 s = GetWorldToScreen(toRl(viewPos), cam3d_);
    return { s.x, s.y, true };
}

int RaylibBackend::ScreenWidth() const  { return GetScreenWidth(); }
int RaylibBackend::ScreenHeight() const { return GetScreenHeight(); }

void RaylibBackend::DrawModel(MeshHandle h, const RMat4& model, const Material& mat) {
    if (h == 0 || h > meshes_.size()) return;
    const wire::GpuMesh& m = meshes_[h - 1];
    const bool additive = mat.blend == BlendMode::Additive;

    // Opaque meshes hide what's behind them; see-through ones are just edges.
    if (!additive && mat.depth_write) wire_.Fill(m, model);

    wire::Shading s;
    s.tint = theme::Remap(mat.color);
    if (mat.lit) { s.heatDir = heatDir_; s.heat = heat_; }
    rlSetLineWidth(additive ? 1.0f : 1.5f);
    if (additive) rlSetBlendMode(RL_BLEND_ADDITIVE);
    if (!mat.depth_write) rlDisableDepthMask();
    wire_.Edges(m, model, s);
    if (!mat.depth_write) rlEnableDepthMask();
    if (additive) rlSetBlendMode(RL_BLEND_ALPHA);
}

void RaylibBackend::DrawLines(const LineVertex* v, size_t count, float width) {
    // Into the display palette. Runs of one colour (a whole path) look it up once.
    remapped_.resize(count);
    RColor in = { 0, 0, 0, 0 }, out = in;
    for (size_t i = 0; i < count; ++i) {
        RColor c = v[i].color;
        if (c.r != in.r || c.g != in.g || c.b != in.b || c.a != in.a) { in = c; out = theme::Remap(c); }
        remapped_[i] = { v[i].pos, out };
    }
    rlSetLineWidth(width > 1.5f ? 1.5f : width);   // thin, crisp strokes
    wire_.Lines(remapped_.data(), count);
}

void RaylibBackend::DrawRocket(const RocketFrame& f) {
    rocket_.Ensure(*this, f.dims);   // build on first sight of each stage config
    heatDir_ = f.vel_dir;
    heat_    = f.heating > 0.03f ? f.heating : 0.0f;
    rocket_.Draw(*this, f);
    heat_    = 0.0f;
}

void RaylibBackend::DrawEarth(const EarthFrame& f) {
    rlSetLineWidth(1.0f);
    int h = GetScreenHeight();
    earth_.Draw(wire_, f, cam_, h > 0 ? (float)GetScreenWidth() / (float)h : 1.0f);
    rlSetLineWidth(1.5f);
    coast_.Draw(wire_, f);
    rlSetLineWidth(1.0f);
    ground_.Draw(*this, f, cam_, labels_);
}

bool RaylibBackend::DrawHud(const HudFrame& hud) {
    HudView v;
    v.width       = GetScreenWidth();
    v.height      = GetScreenHeight();
    v.now         = GetTime();
    v.fps         = GetFPS();
    v.haveHeading = ground_.HaveHeading();
    v.headingDeg  = ground_.HeadingDeg();
    v.labels      = &labels_;
    v.project     = [this](const RVec3& p) { return WorldToScreen(p); };
    hud_.Draw(hud, v);
    return true;
}

void RaylibBackend::DrawRect(int x, int y, int w, int h, RColor c) {
    DrawRectangle(x, y, w, h, toRl(c));
}

void RaylibBackend::DrawRectLines(int x, int y, int w, int h, RColor c) {
    DrawRectangleLines(x, y, w, h, toRl(c));
}

void RaylibBackend::DrawText(const char* text, int x, int y, int font_size, RColor c) {
    hud_.Text(text, x, y, font_size, theme::Remap(c));
}

void RaylibBackend::DrawFPS(int x, int y) { ::DrawFPS(x, y); }

}
