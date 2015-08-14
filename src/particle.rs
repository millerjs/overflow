extern crate rand;

use rand::distributions::{IndependentSample,Range};

use linalg::Matrix;
use linalg::Vec3d;
use overflow::Particle;
use domain::{BOX_X,BOX_Y,BOX_Z,F_G};


#[allow(dead_code)]
impl Particle {

    /// Create a new particle
    pub fn new(rad: f64,
           rx: f64, ry: f64, rz: f64,
           vx: f64, vy: f64, vz: f64) -> Particle{
        let mut rng = rand::thread_rng();
        Particle {
            id: Range::new(0, u64::max_value()).ind_sample(&mut rng),
            r: Vec3d::new(rx, ry, rz),
            v: Vec3d::new(vx, vy, vz),
            a: Vec3d::zero(),
            ap: Vec3d::zero(),
            m: 1.0,
            rad: rad,
            elast: 0.3,
        }
    }

    /// Generate particle with random position
    pub fn random_pos(rad: f64) -> Particle {
        let mut rng = rand::thread_rng();
        let x = Range::new(-BOX_X, BOX_X).ind_sample(&mut rng);
        let y = Range::new(0.0,    BOX_Y).ind_sample(&mut rng);
        let z = Range::new(-BOX_Z, BOX_Z).ind_sample(&mut rng);
        Particle::new(rad, x, y, z, 0.0, 0.0, 0.0)
    }

    /// Get the distance of particle from origin
    pub fn dist(&self) -> f64 {
        (self.r[0]*self.r[0] + self.r[1]*self.r[1] + self.r[2]*self.r[2]).sqrt()
    }

    /// Update positions and velocities using velocity verlet scheme
    pub fn verlet(&mut self, dt: f64) {
        for j in 0..3 {
            self.v[j] += 0.5*(self.a[j] + self.ap[j])*dt;
            self.r[j] += self.v[j]*dt + 0.5*self.a[j]*dt*dt;
        }
    }

    /// Add the acceleration from gravity to the particle
    pub fn force_gravity (&mut self) {
        self.a[1] -= F_G*self.m;
    }

    pub fn surrounding_collisions(&mut self, positions: &Vec<Particle>) {
        for p in positions.iter() {
            let r = self.r - p.r;
            let norm = r.norm();
            let unit_r = [r[0]/norm, r[1]/norm, r[2]/norm];
            let overlap = self.rad + p.rad - norm;
            if norm > 0.0 && overlap > 0.0 {
                for i in 0..3 {
                    self.a[i] += 5e3 * overlap * unit_r[i];
                    self.a[i] -= 1.0 * self.v[i];
                }
            }
        }
    }

    /// Call force functions to update the acceleration on the
    /// particle at this point in time
    pub fn force_total (&mut self, positions: &Vec<Particle>) {
        self.ap = self.a;
        self.a = Vec3d::zero();
        self.surrounding_collisions(positions);
        self.force_gravity();
    }

    pub fn wall_collision (&mut self, dir: usize) {
        self.v[dir] = - self.elast * self.v[dir];
        self.v[(dir + 1) % 3] = (1.0-1e-4) * self.v[(dir + 1) % 3];
        self.v[(dir + 2) % 3] = (1.0-1e-4) * self.v[(dir + 2) % 3];
    }

    /// Reflect particles off bounding box
    pub fn check_boundaries (&mut self) {
        let dim = [BOX_X, BOX_Y, BOX_Z];
        // X
        if self.r[0] < -dim[0]  { self.r[0] = -dim[0]; self.wall_collision(0); }
        if self.r[0] >  dim[0]  { self.r[0] =  dim[0]; self.wall_collision(0); }

        // Y
        let y_min = - (BOX_Y as f64 / BOX_X as f64 / 2.0) *
            self.r[0] - BOX_Y as f64 / 2.0;
        let vp = &Matrix::householder(
            Vec3d::new(1.0, 1.0, 0.0), 1.0-1e-4)*self.v;
        if self.r[1] <  y_min   { self.r[1] =  y_min;  self.v = vp; }
        if self.r[1] >  dim[1]  { self.r[1] =  dim[1]; self.wall_collision(1); }

        // Z
        if self.r[2] < -dim[2]  { self.r[2] = -dim[2]; self.wall_collision(2); }
        if self.r[2] >  dim[2]  { self.r[2] =  dim[2]; self.wall_collision(2); }
    }

    /// Perform all updates on particle given neighbors and delta t
    pub fn update(&mut self, positions: &Vec<Particle>, dt: f64) {
        self.force_total(positions);
        self.verlet(dt);
        self.check_boundaries();
    }

}
