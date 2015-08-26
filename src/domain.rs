extern crate rand;

use threadpool::ThreadPool;
use rand::distributions::{IndependentSample,Range};
use std::sync::mpsc::channel;
use std::sync::Arc;

use linalg::Vec3d;
use overflow::{Domain,Particle};

// Bounding box dim constatns
pub const BOX_X: f64 = 600.0;
pub const BOX_Y: f64 = 600.0;
pub const BOX_Z: f64 = 300.0;
pub const BOX: [f64; 3] = [BOX_X, BOX_Y, BOX_Z];
pub const F_G: f64 = 1.0e3;

#[allow(dead_code)]
impl Domain {

    pub fn add_particles_grid(&mut self, rad: f64, n: u32) {
        let mut rng = rand::thread_rng();
        let range = Range::new(-rad/5.0, rad/5.0);
        let dist = 2.0 * rad;
        for x in 1..n {
            for y in 1..n {
                for z in 1..5 {
                    self.particles.push(Particle::new(
                        rad,
                        x as f64 * dist + range.ind_sample(&mut rng),
                        y as f64 * dist + range.ind_sample(&mut rng),
                        z as f64 * dist + range.ind_sample(&mut rng),
                        0.0, 0.0, 0.0));
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn add_particles_random(&mut self, rad: f64, count: usize) {
        for _ in 0..count {
            self.particles.push(Particle::random_pos(rad));
        }
    }

    pub fn draw_bounding_box(&mut self) {
        let p = [1.0; 4];
        let lft_y = -BOX_Y;
        // top
        self.camera.draw_line(Vec3d::new( BOX_X,  BOX_Y,  BOX_Z), Vec3d::new( BOX_X,  BOX_Y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new(-BOX_X,  BOX_Y,  BOX_Z), Vec3d::new(-BOX_X,  BOX_Y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new( BOX_X,  BOX_Y, -BOX_Z), Vec3d::new(-BOX_X,  BOX_Y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new( BOX_X,  BOX_Y,  BOX_Z), Vec3d::new(-BOX_X,  BOX_Y,  BOX_Z), p);
        // bottom
        self.camera.draw_line(Vec3d::new( BOX_X, -BOX_Y,  BOX_Z), Vec3d::new( BOX_X, -BOX_Y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new(-BOX_X,  lft_y,  BOX_Z), Vec3d::new(-BOX_X,  lft_y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new( BOX_X, -BOX_Y, -BOX_Z), Vec3d::new(-BOX_X,  lft_y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new( BOX_X, -BOX_Y,  BOX_Z), Vec3d::new(-BOX_X,  lft_y,  BOX_Z), p);
        // sides
        self.camera.draw_line(Vec3d::new( BOX_X, -BOX_Y,  BOX_Z), Vec3d::new( BOX_X,  BOX_Y,  BOX_Z), p);
        self.camera.draw_line(Vec3d::new(-BOX_X,  lft_y,  BOX_Z), Vec3d::new(-BOX_X,  BOX_Y,  BOX_Z), p);
        self.camera.draw_line(Vec3d::new( BOX_X, -BOX_Y, -BOX_Z), Vec3d::new( BOX_X,  BOX_Y, -BOX_Z), p);
        self.camera.draw_line(Vec3d::new(-BOX_X,  lft_y, -BOX_Z), Vec3d::new(-BOX_X,  BOX_Y, -BOX_Z), p);
    }

    pub fn render(&mut self) {
        self.camera.clear([0.0; 4]);
        self.draw_bounding_box();
        for p in self.particles.iter(){
            self.camera.draw_particle(p, [0.5; 4]);
        }
        self.camera.save();
    }

    pub fn update(&mut self) {
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

    pub fn print_state(&self) {
        println!("\tDomain:\t\trender_step={}\t\tt={:.*}",
                 self.camera.render_step, 2, self.t);
    }

}
