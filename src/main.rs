extern crate num;
extern crate rand;
extern crate image;

use std::f64::consts::PI;
use std::thread;
use rand::distributions::{IndependentSample, Range};
use image::{ImageBuffer, Rgb};

const FOVY: f64 = 60.0;

const BOX_X: f64 = 600.0;
const BOX_Y: f64 = 300.0;
const BOX_Z: f64 = 300.0;

const IMAGE_X: i32 = 400;
const IMAGE_Y: i32 = 300;

const SCREEN: f64 = 600.0;
const SCALE: f64 = IMAGE_X as f64 / 1200.0;

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
    width: i32,
    height: i32,
}

#[allow(dead_code)]
impl Image {
    fn new(width: i32, height: i32) -> Image {
        Image {
            pixels: vec![[0; 4]; width as usize * height as usize],
            width: width,
            height: height,
        }
    }

    fn set_pixel(&mut self, x: i32, y: i32, pixel: [f32; 4]) {
        // println!("{} {}", x, y);
        if x < self.width && x > 0 && y < self.height && y > 0 {
            for i in 0..4 {
                let idx = (self.height * x + y) as usize;
                self.pixels[idx][i] = (pixel[i]*255.0) as u8;
            }
        }
    }

    fn get_pixel(&mut self, x: i32, y: i32) -> [u8; 4]{
        self.pixels[(self.height * x + y) as usize]
    }

    fn set_point(&mut self, x: f64, y: f64, pixel: [f32; 4]) {
        let x = self.translate(x);
        let y = self.translate(y);
        self.set_pixel(x, y, pixel);
    }

    fn translate(&self, v: f64) -> i32 {
         (v * SCALE) as i32
    }

    fn draw_circle(&mut self, _x: f64, _y: f64, r: f64, pixel: [f32; 4]) {
        let x0 = (self.width/2 - self.translate(_x)) as i32;
        let y0 = (self.height/2 - self.translate(_y)) as i32;
        let mut x = r as i32;
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

    fn clear(&mut self, color: [f32; 4]) {
        for i in 0..(self.width*self.height) as usize {
            for j in 0..4 {
                self.pixels[i][j] = (color[j]*255 as f32) as u8;
            }
        }
    }

    fn save(&mut self, step: u32) -> thread::JoinHandle<()> {
        println!(" - Generating image at step {}", step);
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
            println!(" - Saving image to '{}'", path);
            imgbuf.save(&*path).unwrap();
            println!(" - Wrote '{}'", path);
        })

    }

}


pub struct Domain {
    render_step: u32,
    t: f64,
    dt: f64,
    particles: Vec<Particle>,
    bbox: Vec<Particle>,
    image: Image,
}

#[allow(dead_code)]
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

#[allow(dead_code)]
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
        (1.0 / fov.tan() * self.rad /
            (d*d - self.rad*self.rad).sqrt() * (IMAGE_Y/2) as f64).abs()
    }

    /// Get the x or y location as projected on the screen
    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/(self.r[2] + SCREEN*2.0)
    }

    /// Update positions and velocities using velocity verlet scheme
    fn verlet(&mut self, dt: f64){
        for j in 0..3 {
            self.r[j] += self.v[j]*dt + 0.5*self.a[j]*dt*dt;
            self.v[j] += 0.5*(self.a[j] + self.ap[j])*dt;
            self.v[j] = self.v[j].min(1e4).max(-1e4);
        }
    }

    /// Add the acceleration from gravity to the particle
    fn force_gravity (&mut self) {
        self.a[1] -= 1000.0*self.m;
    }

    /// Add the acceleration from pseudo Lennard-Jones potentials from
    /// surrounding particles
    fn force_lj (&mut self, positions: &Vec<[f64; 3]>) {
        let sigma = self.rad / 1.122 * 2.0;
        let eps = 1.0e5;
        for p in positions.iter() {
            let r = v_sub(self.r, *p);
            let r_sq = v_dot(r, r);
            if r_sq > 1e-5 {
                let norm = r_sq.sqrt();
                let temp_7 = (sigma / norm).powi(7);
                let temp_13 = (sigma / norm).powi(13);
                let a = eps / sigma * (2.0 * temp_13 - temp_7);
                for j in 0..3 {
                    self.a[j] -= a * r[0] / norm;
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

        let n = 5;
        let rad = 10.0;
        let mut rng = rand::thread_rng();
        let range = Range::new(-rad/3.0, rad/3.0);

        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    self.particles.push(Particle::new(
                        rad,
                        x as f64 * rad * 5.0 + range.ind_sample(&mut rng),
                        y as f64 * rad * 5.0 + range.ind_sample(&mut rng),
                        z as f64 * rad * 5.0 + range.ind_sample(&mut rng),
                        0.0, 0.0, 0.0));
                }
            }
        }

        // for _ in 0..250 {
        //     self.particles.push(Particle::random_pos());
        // }

        // self.particles.push(Particle::new(
        //     10.0, -BOX_X/2.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        // self.particles.push(Particle::new(
        //     10.0, BOX_X/2.0, 0.0,  0.0, 0.0, 0.0, 0.0));


        self.bbox.push(Particle::new(10.0, -BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0, -BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0, -BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0,  BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
        self.bbox.push(Particle::new(10.0, -BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
    }


    fn render(&mut self) {
        self.image.clear([0.0; 4]);
        for p in self.bbox.iter(){
            self.image.draw_circle(
                p.projected(p.r[0]), p.projected(p.r[1]), p.projected_radius(),
                [0.5, 0.5, 0.8, 1.0]);
        }
        for p in self.particles.iter(){
            self.image.draw_circle(
                p.projected(p.r[0]), p.projected(p.r[1]), p.projected_radius(),
                [0.5; 4]);

        }
        self.render_step += 1;
    }

    fn update(&mut self) {
        let mut positions: Vec<[f64; 3]> = vec![];
        for p in self.particles.iter(){
            positions.push(p.r);
        }
        for p in self.particles.iter_mut(){
            p.update(&positions, self.dt);
        }
        self.t += self.dt;
    }

    fn print_state(&self) {
        println!("Domain:\tt={:.*}", 2, self.t);
    }

}


fn main() {

    let mut domain = Domain {
        particles: vec![],
        bbox: vec![],
        render_step: 0,
        dt: 1.0e-2,
        t: 0.0,
        image: Image::new(IMAGE_X, IMAGE_Y),
    };

    println!("Created image buffer of length {} KB", IMAGE_X*IMAGE_Y/1000);
    domain.setup();

    let max_time = 3.0;
    let mut step = 0;

    while domain.t < max_time {
        domain.update();
        if step % (1.0/domain.dt/10.0) as i32 == 0 {
            domain.render();
            domain.image.save(domain.render_step);
            domain.print_state();
        }
        step += 1;
    }

    thread::sleep_ms(1000);

}
