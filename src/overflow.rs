use linalg::Vec3d;

#[derive(Debug, Copy, Clone)]
pub struct Particle {
    pub id: u64,
    pub r: Vec3d,      // location
    pub v: Vec3d,      // velocity
    pub a: Vec3d,      // acceleration
    pub ap: Vec3d,     // acceleration from previous time-step
    pub m: f64,        // mass
    pub rad: f64,      // radius
    pub elast: f64,
    pub dens: f64,
}


pub struct Camera {
    pub render_step: u32,
    pub pixels: Vec<[u8; 4]>,
    pub screen: f64,
    pub fovd: f64,
    pub width: i32,
    pub height: i32,
    pub r: Vec3d,
    pub theta: Vec3d,  // camera angle
}


pub struct Domain {
    pub t: f64,
    pub dt: f64,
    pub particles: Vec<Particle>,
    pub bbox: Vec<Particle>,
    pub camera: Camera,
}
