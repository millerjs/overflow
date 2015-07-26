extern crate num;
extern crate rand;
extern crate image;

use std::fs::File;
use std::path::Path;
use std::f64::consts::PI;
use std::cmp;
use rand::Rng;
use rand::distributions::{IndependentSample, Range};
use image::ImageBuffer;

const BACKGROUND: [f32; 4] = [0.0, 0.0, 0.0, 1.0];
const FOVY: f64 = 60.0;
const BOX_X: f64 = 600.0;
const BOX_Y: f64 = 300.0;
const BOX_Z: f64 = 300.0;
const IMAGE_X: usize = 1200;
const IMAGE_Y: usize = 800;
const SCREEN: f64 = 600.0;
const WINDOW: [i64; 2] = [1200, 800];
const SCALE: f64 = 1.0;

pub struct Particle {
    r: [f64; 3],      // location
    v: [f64; 3],      // velocity
    a: [f64; 3],      // acceleration
    ap: [f64; 3],     // acceleration from previous time-step
    m: f64,           // mass
    rad: f64,         // radius
    elast: f64,
}

pub struct Image {
    pixels: Vec<[u8; 4]>,
    width: usize,
    height: usize,
}

impl Image {
    fn new(width: usize, height: usize) -> Image {
        Image {
            pixels: vec![[0; 4]; width as usize * height as usize],
            width: width,
            height: height,
        }
    }

    fn set_pixel(&mut self, x: usize, y: usize, pixel: [f32; 4]) {
        if x < self.width && x > 0 && y < self.height && y > 0 {
            for i in 0..5 {
                self.pixels[self.width * x + y][i] = (pixel[i]*255.0) as u8;
            }
        }
    }

    fn set_point(&mut self, x: f64, y: f64, pixel: [f32; 4]) {
        let t = self.translate(x, y);
        self.set_pixel(t[0], t[1], pixel);
    }

    fn translate(&self, x: f64, y: f64) -> [usize; 2] {
        [self.width / 2 - (x * SCALE) as usize,
         self.height / 2 - (x * SCALE) as usize]
    }

    fn draw_circle(&mut self, x: f64, y: f64, r: f64, pixel: [f32; 4]) {
        let t = self.translate(x, y);
        let x = t[0] as u32;
        let y = t[1] as u32;
        let x0 = x;
        while x < x0 {
            let x = ((x*x - 2*y - 1) as f64).sqrt();
            let y = ((r*r - x*x) as f64).sqrt();
            self.set_pixel(x as usize, y as usize, pixel);
        }
    }

    fn clear(&mut self, color: [f32; 4]) {
        for i in 0..self.width*self.height {
            self.pixels[i] = [0; 4];
        }
    }

    fn save(mut self, path: &str) {
        let mut imgbuf = image::ImageBuffer::new(
            self.width as u32, self.height as u32);

        // Iterate over the coordiantes and pixels of the image
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            // let cy = y as f32 * scaley - 2.0;
            // let cx = x as f32 * scalex - 2.0;

            // let mut z = Complex::new(cx, cy);
            // let c = Complex::new(-0.4, 0.6);

            // let mut i = 0;

            // for t in (0..max_iterations) {
            //     if z.norm() > 2.0 {
            //         break
            //     }
            //     z = z * z + c;
            //     i = t;
            // }

            // // Create an 8bit pixel of type Luma and value i
            // // and assign in to the pixel at position (x, y)
            // *pixel = image::Luma([i as u8]);
        }


        let ref mut fout = File::create(&Path::new(path)).unwrap();
        let _ = image::ImageLuma8(imgbuf).save(fout, image::PNG);
    }

}


pub struct Domain {
    step: i32,
    t: f64,
    dt: f64,
    particles: Vec<Particle>,
    bbox: Vec<Particle>,
    image: Image,
}

fn v_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for i in 0..3 {
        c[i] = a[i] + b[i];
    };
    c
}

fn v_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for i in 0..3 {
        c[i] = a[i] - b[i];
    };
    c
}

fn v_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut c = 0.0;
    for i in 0..3 {
        c += a[i] * b[i];
    };
    c
}


impl Particle {

    /// Create a new particle
    fn new(rad: f64,
           rx: f64, ry: f64, rz: f64,
           vx: f64, vy: f64, vz: f64) -> Particle{
        Particle {
            r: [rx, ry, rz],
            v: [vx, vy, vz],
            a: [0.0, 0.0, 0.0],
            ap: [0.0, 0.0, 0.0],
            m: 1.0,
            rad: rad,
            elast: 0.4,
        }
    }

    /// Generate particle with random position
    fn random_pos() -> Particle {
        let mut rng = rand::thread_rng();
        let x = Range::new(-BOX_X, BOX_X).ind_sample(&mut rng);
        let y = Range::new(-BOX_Y, BOX_Y).ind_sample(&mut rng);
        let z = Range::new(-BOX_Z, BOX_Z).ind_sample(&mut rng);
        Particle::new(10.0, x, y, z, 100.0, 0.0, 0.0)
    }

    /// Get the distance of particle from origin
    fn dist(&self) -> f64 {
        (self.r[0]*self.r[0] + self.r[1]*self.r[1] + self.r[2]*self.r[2]).sqrt()
    }

    /// Scale the apparent radius of the ellipse for perspective
    fn projected_radius(&self) -> f64 {
        let fov = FOVY / 2.0 * PI / 180.0;
        let d = self.r[2] + SCREEN;
        1.0 / fov.tan() * self.rad /
            (d*d - self.rad*self.rad).sqrt() * WINDOW[0] as f64/1200.0
    }

    /// Get the x or y location as projected on the screen
    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/(self.r[2] + SCREEN*2.0) * WINDOW[0] as f64/1200.0
    }

    /// Update positions and velocities using velocity verlet scheme
    fn verlet(&mut self, dt: f64){
        for j in 0..3 {
            self.r[j] += self.v[j]*dt + 0.5*self.a[j]*dt*dt;
            self.v[j] += 0.5*(self.a[j] + self.ap[j])*dt;
        }
    }

    /// Add the acceleration from gravity to the particle
    fn force_gravity (&mut self) {
        self.a[1] -= 1000.0*self.m;
    }

    /// Add the acceleration from pseudo Lennard-Jones potentials from
    /// surrounding particles
    fn force_lj (&mut self, positions: &Vec<[f64; 3]>) {
        let sigma = self.rad * 1.0e-5;
        let eps = 100000.0;
        for p in positions.iter() {
            let r = v_sub(self.r, *p);
            let r_sq = v_dot(r, r);
            if r_sq.abs() > 0.01 && r_sq.abs() < 5.0 {
                for j in 0..3 {
                    let temp_6 = (sigma / self.r[j]).powi(6);
                    let temp_12 = temp_6 * temp_6;
                    self.a[j] = (eps * r[j] * (temp_12 - 0.5 * temp_6)/
                                  (r_sq * self.m));
                }
            }
        }
    }

    /// Call force functions to update the acceleration on the
    /// particle at this point in time
    fn force_total (&mut self, positions: &Vec<[f64; 3]>) {
        self.ap = self.a;
        self.a = [0.0; 3];
        self.force_lj(positions);
        self.force_gravity();
    }

    fn wall_collision (&mut self) {
        for i in 0..3 {
            self.v[i] = self.elast * self.v[i];
        }
    }

    /// Reflect particles off bounding box
    fn check_boundaries (&mut self) {
        let buffer = 0.0;
        if self.r[0] < -BOX_X  {
            self.v[0] = - self.v[0];
            self.r[0] = -BOX_X + buffer;
            self.wall_collision();
        }
        if self.r[0] > BOX_X  {
            self.v[0] = -self.v[0];
            self.r[0] = BOX_X - buffer;
            self.wall_collision();
        }
        if self.r[1] < -BOX_Y  {
            self.v[1] = -self.v[1];
            self.r[1] = -BOX_Y + buffer;
            self.wall_collision();
        }
        if self.r[1] > BOX_Y  {
            self.v[1] = -self.v[1];
            self.r[1] = BOX_Y -  buffer;
            self.wall_collision();
        }
        if self.r[2] < -BOX_Z  {
            self.v[2] = -self.v[2];
            self.r[2] = -BOX_Z + buffer;
            self.wall_collision();
        }
        if self.r[2] > BOX_Z  {
            self.v[2] = -self.v[2];
            self.r[2] = BOX_Z - buffer;
            self.wall_collision();
        }
    }

    /// Perform all updates on particle given neighbors and delta t
    fn update(&mut self, positions: &Vec<[f64; 3]>, dt: f64) {
        self.force_total(positions);
        self.verlet(dt);
        self.check_boundaries();
    }

}


impl Domain {

    fn setup(&mut self) {

        for i in 0..100 {
            self.particles.push(Particle::random_pos());
        }

        self.bbox.push(Particle::new(20.0, -BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0,  BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0,  BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0, -BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0, -BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0,  BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0,  BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(20.0, -BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
    }

    fn render(&mut self) {

    }

    fn update(&mut self) {
        let mut positions: Vec<[f64; 3]> = vec![];
        for p in self.particles.iter(){
            positions.push(p.r);
        }
        for p in self.particles.iter_mut(){
            p.update(&positions, self.dt);
        }
    }

}


fn main() {

    let max_iterations = 256u16;

    let mut domain = Domain {
        particles: vec![],
        bbox: vec![],
        step: 0,
        dt: 1.0,
        t: 0.0,
        image: Image::new(IMAGE_X, IMAGE_Y),
    };

    println!("Created image buffer of length {} KB", IMAGE_X*IMAGE_Y/1000);

    domain.setup();

    for i in 0..5 {
        domain.update();
    }

    domain.image.draw_circle(0.0, 0.0, 20.0, [0.5; 4]);
    domain.image.save("output/test.png");

}
