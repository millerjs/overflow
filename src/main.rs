extern crate rand;
extern crate image;
extern crate threadpool;

// Usage statements
use std::f64::consts::PI;
use std::thread;
use rand::distributions::{IndependentSample, Range};
use image::{ImageBuffer, Rgb};
use threadpool::ThreadPool;
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::cmp::{min, max};

mod linalg;
use linalg::Matrix;
use linalg::Vec3d;

// Bounding box dim constatns
const BOX_X: f64 = 300.0;
const BOX_Y: f64 = 600.0;
const BOX_Z: f64 = 200.0;

const F_G: f64 = 1000.0;

#[derive(Debug, Copy, Clone)]
pub struct Particle {
    id: u64,
    r: Vec3d,      // location
    v: Vec3d,      // velocity
    a: Vec3d,      // acceleration
    ap: Vec3d,     // acceleration from previous time-step
    m: f64,        // mass
    rad: f64,      // radius
    elast: f64,
}

pub struct Camera {
    pixels: Vec<[u8; 4]>,
    screen: f64,
    fovd: f64,
    width: i32,
    height: i32,
    camera: Vec3d,
    theta: Vec3d,  // camera angle
}

pub struct Domain {
    render_step: u32,
    t: f64,
    dt: f64,
    particles: Vec<Particle>,
    bbox: Vec<Particle>,
    camera: Camera,
}

#[allow(dead_code)]
impl Camera {

    /// Scale the apparent radius of the ellipse for perspective
    fn projected_radius(&self, r: Vec3d, rad: f64) -> f64 {
        let fov = self.fovd / 2.0 * PI / 180.0;
        let d = (r - self.camera).norm();
        (1.2e-1 / fov.tan() * rad /
            (d*d - rad*rad).sqrt() * (self.width) as f64).abs()
    }

    /// Get the x, y location as projected on the screen
    fn projected(&self, v: Vec3d) -> [i32; 2] {
        let a = Matrix::rot_x(-self.theta[0]);
        let b = Matrix::rot_y(-self.theta[1]);
        let c = Matrix::rot_z(-self.theta[2]);
        let d = &(&(&a * &b) * &c) * (v-self.camera);
        let bx = (d[0] * self.width  as f64)/(d[2]*self.screen);
        let by = (d[1] * self.height as f64)/(d[2]*self.screen);
        [self.width/2 - bx as i32, self.height/2 - by as i32]
    }

    fn new(width: i32, height: i32, camera: Vec3d) -> Camera {
        Camera {
            pixels: vec![[0; 4]; width as usize * height as usize],
            width: width,
            height: height,
            camera: camera,
            screen: 5.0,
            fovd: 60.0,
            theta: Vec3d::new(0.0, 0.0, 0.0),
        }
    }

    fn set_pixel(&mut self, x: i32, y: i32, pixel: [f32; 4]) {
        if x < self.width && x > 0 && y < self.height && y > 0 {
            for i in 0..4 {
                let idx = (self.height * x + y) as usize;
                self.pixels[idx][i] = (pixel[i]*255.0) as u8;
            }
        }
    }

    fn get_pixel(&mut self, x: i32, y: i32) -> [u8; 4]{
        if x < self.width && x > 0 && y < self.height && y > 0 {
            let idx = (self.height * x + y) as usize;
            self.pixels[idx]
        } else {
            [0; 4]
        }
    }

    fn draw_particle(&mut self, p: &Particle, pixel: [f32; 4]) {
        let translated = self.projected(p.r);
        let (x0, y0) = (translated[0], translated[1]);
        let mut x = self.projected_radius(p.r, p.rad) as i32;
        let mut y = 0;
        let mut condition = 1 - x;

        while x > y {
            self.set_pixel(( x + x0) as i32, ( y + y0) as i32, pixel);
            self.set_pixel(( y + x0) as i32, ( x + y0) as i32, pixel);
            self.set_pixel((-x + x0) as i32, ( y + y0) as i32, pixel);
            self.set_pixel((-y + x0) as i32, ( x + y0) as i32, pixel);
            self.set_pixel((-x + x0) as i32, (-y + y0) as i32, pixel);
            self.set_pixel((-y + x0) as i32, (-x + y0) as i32, pixel);
            self.set_pixel(( x + x0) as i32, (-y + y0) as i32, pixel);
            self.set_pixel(( y + x0) as i32, (-x + y0) as i32, pixel);
            y += 1;
            if condition <=0 {
                condition += 2 * y + 1;
            } else {
                x -= 1;
                condition += 2 * (y - x) + 1;
            }
        }
    }

    fn draw_line(&mut self, _p0: Vec3d, _p1: Vec3d, pixel: [f32; 4]) {
        let p0 = self.projected(_p0);
        let p1 = self.projected(_p1);
        let (x0, y0, x1, y1) = (p0[0], p0[1], p1[0], p1[1]);
        let (dx, dy) = (x1 - x0, y1 - y0);
        if dx < dy {
            for x in (min(x0, x1)..max(x0, x1)) {
                self.set_pixel(x, y0+(x-x0)*dy/dx, pixel);
            }
        } else {
            for y in (min(y0, y1)..max(y0, y1)) {
                self.set_pixel(x0+(y-y0)*dx/dy, y, pixel);
            }
        }
    }


    fn clear(&mut self, color: [f32; 4]) {
        for i in 0..(self.width*self.height) as usize {
            for j in 0..4 {
                self.pixels[i][j] = (color[j]*255 as f32) as u8;
            }
        }
    }

    fn save(&mut self, step: u32) -> thread::JoinHandle<()> {
        let mut imgbuf = ImageBuffer::<Rgb<u8>>::new(
            self.width as u32, self.height as u32);

        // Iterate over the coordiantes and pixels of the image
        for x in 0..self.width {
            for y in 0..self.height {
                let px = self.get_pixel(x, y);
                imgbuf.get_pixel_mut(x as u32, y as u32).data = [
                    px[0], px[1], px[2]];
            }
        }

        thread::spawn(move || {
            let path = format!("output/test_{:08}.png", step);
            imgbuf.save(&*path).unwrap();
        })

    }
}


#[allow(dead_code)]
impl Particle {

    /// Create a new particle
    fn new(rad: f64,
           rx: f64, ry: f64, rz: f64,
           vx: f64, vy: f64, vz: f64) -> Particle{
        let mut rng = rand::thread_rng();
        Particle {
            id: Range::new(0, std::u64::MAX).ind_sample(&mut rng),
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
    fn random_pos() -> Particle {
        let mut rng = rand::thread_rng();
        let x = Range::new(-BOX_X, 0.0).ind_sample(&mut rng);
        let y = Range::new(-BOX_Y, 0.0).ind_sample(&mut rng);
        let z = Range::new(-BOX_Z, 0.0).ind_sample(&mut rng);
        Particle::new(10.0, x, y, z, 0.0, 0.0, 0.0)
    }

    /// Get the distance of particle from origin
    fn dist(&self) -> f64 {
        (self.r[0]*self.r[0] + self.r[1]*self.r[1] + self.r[2]*self.r[2]).sqrt()
    }

    /// Update positions and velocities using velocity verlet scheme
    fn verlet(&mut self, dt: f64) {
        for j in 0..3 {
            self.v[j] += 0.5*(self.a[j] + self.ap[j])*dt;
            self.r[j] += self.v[j]*dt + 0.5*self.a[j]*dt*dt;
        }
    }

    /// Add the acceleration from gravity to the particle
    fn force_gravity (&mut self) {
        self.a[1] -= F_G*self.m;
    }

    fn surrounding_collisions(&mut self, positions: &Vec<Particle>) {
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
    fn force_total (&mut self, positions: &Vec<Particle>) {
        self.ap = self.a;
        self.a = Vec3d::zero();
        self.surrounding_collisions(positions);
        self.force_gravity();
    }

    fn wall_collision (&mut self, dir: usize) {
        self.v[dir] = - self.elast * self.v[dir];
        self.v[(dir + 1) % 3] = (1.0-1e-4) * self.v[(dir + 1) % 3];
        self.v[(dir + 2) % 3] = (1.0-1e-4) * self.v[(dir + 2) % 3];
    }

    /// Reflect particles off bounding box
    fn check_boundaries (&mut self) {
        let dim = [BOX_X, BOX_Y, BOX_Z];
        for i in 0..3 {
            if self.r[i] < -dim[i]  {
                self.r[i] = -dim[i];
                self.wall_collision(i);
            }
            if self.r[i] > dim[i]  {
                self.r[i] = dim[i];
                self.wall_collision(i);
            }
        }
    }

    /// Perform all updates on particle given neighbors and delta t
    fn update(&mut self, positions: &Vec<Particle>, dt: f64) {
        self.force_total(positions);
        self.verlet(dt);
        self.check_boundaries();
    }

}


impl Domain {

    #[allow(dead_code)]
    fn add_particles_grid(&mut self, rad: f64, n: u32) {
        let mut rng = rand::thread_rng();
        let range = Range::new(-rad/5.0, rad/5.0);
        for x in 1..n {
            for y in 1..n {
                for z in 1..5 {
                    self.particles.push(Particle::new(
                        rad,
                        - BOX_X + (
                            x as f64 * rad * 2.1 + range.ind_sample(&mut rng)),
                        - (y as f64 * rad * 2.1 + range.ind_sample(&mut rng)),
                        - BOX_Z + (
                            z as f64 * rad * 2.1 + range.ind_sample(&mut rng)),
                        0.0, 0.0, 0.0));
                }
            }
        }
    }

    #[allow(dead_code)]
    fn add_particles_random(&mut self, count: usize) {
        for _ in 0..count {
            self.particles.push(Particle::random_pos());
        }
    }

    fn add_bounding_box(&mut self) {
        self.bbox.push(Particle::new(10.0, -BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0, -BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0, -BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0, -BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
    }

    fn setup(&mut self) {
        self.add_particles_grid(10.0, 4);
        self.add_bounding_box();
    }

    fn draw_bounding_box(&mut self) {
        let p = [1.0; 4];
        for z in [BOX_Z, -BOX_Z].iter() {
            self.camera.draw_line(
                Vec3d::new(-BOX_X, BOX_Y, *z),
                Vec3d::new(BOX_X,  BOX_Y, *z), p);
            self.camera.draw_line(
                Vec3d::new(-BOX_X, -BOX_Y, *z),
                Vec3d::new(BOX_X, -BOX_Y, *z), p);
            self.camera.draw_line(
                Vec3d::new(BOX_X, -BOX_Y, *z),
                Vec3d::new(BOX_X, BOX_Y, *z), p);
            self.camera.draw_line(
                Vec3d::new(-BOX_X, -BOX_Y, *z),
                Vec3d::new(-BOX_X, BOX_Y, *z), p);

        }
        for y in [BOX_Y, -BOX_Y].iter() {
            self.camera.draw_line(
                Vec3d::new(-BOX_X, *y, -BOX_Z),
                Vec3d::new(-BOX_X, *y, BOX_Z), p);
            self.camera.draw_line(
                Vec3d::new(BOX_X, *y, -BOX_Z),
                Vec3d::new(BOX_X, *y, BOX_Z), p);
        }
    }

    fn render(&mut self) {
        self.camera.clear([0.0; 4]);
        self.draw_bounding_box();
        for p in self.bbox.iter(){
            self.camera.draw_particle(p, [0.5, 0.5, 0.8, 1.0]);
        }
        for p in self.particles.iter(){
            self.camera.draw_particle(p, [0.5; 4]);
        }
        self.render_step += 1;
    }

    fn update(&mut self) {
        let mut positions: Vec<Particle> = vec![];
        for p in self.particles.iter(){
            positions.push(p.clone());
        }

        let shared_pos = Arc::new(positions);
        let pool = ThreadPool::new(32);

        let (tx, rx) = channel();
        let dt = self.dt;

        for p in self.particles.iter_mut() {
            let tx = tx.clone();
            let pos = shared_pos.clone();
            let mut _p = p.clone();
            pool.execute(move|| {
                _p.update(&pos, dt);
                let _ = tx.send(_p);
            });
        }

        for i in 0..self.particles.len() {
            self.particles[i] = rx.recv().unwrap();
        }

        self.t += self.dt;
    }

    fn print_state(&self) {
        println!("\tDomain:\t\trender_step={}\t\tt={:.*}",
                 self.render_step, 2, self.t);
    }

}


fn main() {

    let mut domain = Domain {
        particles: vec![],
        bbox: vec![],
        render_step: 0,
        dt: 1.0e-2,
        t: 0.0,
        camera: Camera::new(800, 500, Vec3d::new(0.0, 0.0, -600.0)),
    };

    println!("Created image buffer of length {} KB",
             domain.camera.width*domain.camera.height/1000);
    domain.setup();

    let max_time = 1.0;
    let mut step = 0;

    while domain.t < max_time {
        domain.update();
        if step % 5 == 0 {
            domain.print_state();
            domain.render();
            domain.camera.save(domain.render_step);
        }
        step += 1;
    }

    thread::sleep_ms(2000);

}
