extern crate piston;
extern crate graphics;
extern crate glutin_window;
extern crate opengl_graphics;
extern crate rand;

use piston::window::WindowSettings;
use piston::event::*;
use glutin_window::GlutinWindow as Window;
use opengl_graphics::{ GlGraphics, OpenGL };
use graphics::*;
use std::f64::consts::PI;
use std::ops::{Sub, Add};
use rand::Rng;
use rand::distributions::{IndependentSample, Range};

const BACKGROUND: [f32; 4] = [0.0, 0.0, 0.0, 1.0];
const FOVY: f64 = 60.0;
const BOX_X: f64 = 600.0;
const BOX_Y: f64 = 300.0;
const BOX_Z: f64 = 300.0;
const SCREEN: f64 = 600.0;
const WINDOW: [i64; 2] = [1200, 800];

pub struct Particle {
    r: [f64; 3],      // location
    v: [f64; 3],      // velocity
    a: [f64; 3],      // acceleration
    ap: [f64; 3],     // acceleration from previous time-step
    m: f64,           // mass
    rad: f64,         // radius
    elast: f64,
    ellipse: Ellipse, // Graphics
}

pub struct App {
    gl: GlGraphics,
    particles: Vec<Particle>,
    bbox: Vec<Particle>,
}

impl Particle {
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
            ellipse: Ellipse::new([0.2, 0.2, 0.2, 1.0]),
        }
    }

    fn random_pos() -> Particle {
        let mut rng = rand::thread_rng();
        let x = Range::new(-BOX_X, BOX_X).ind_sample(&mut rng);
        let y = Range::new(-BOX_Y, BOX_Y).ind_sample(&mut rng);
        let z = Range::new(-BOX_Z, BOX_Z).ind_sample(&mut rng);
        Particle::new(20.0, x, y, z, 0.0, 0.0, 0.0)
    }

    fn dist(&self) -> f64 {
        (self.r[0]*self.r[0] + self.r[1]*self.r[1] + self.r[2]*self.r[2]).sqrt()
    }

    fn projected_radius(&self) -> f64 {
        let fov = FOVY / 2.0 * PI / 180.0;
        let d = self.r[2] + SCREEN;
        1.0 / fov.tan() * self.rad /
            (d*d - self.rad*self.rad).sqrt() * WINDOW[0] as f64/1200.0
    }

    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/(self.r[2] + SCREEN*2.0) * WINDOW[0] as f64/1200.0
    }

    fn verlet(&mut self, dt: f64){
        for j in 0..3 {
            self.r[j] += self.v[j]*dt + 0.5*self.a[j]*dt*dt;
            self.v[j] += 0.5*(self.a[j] + self.ap[j])*dt;
        }
    }

    fn force_gravity (&mut self) {
        self.a[1] -= 1000.0*self.m;
    }

    fn force_total (&mut self) {
        self.ap = self.a;
        self.a = [0.0; 3];
        self.force_gravity();
    }

    fn check_boundaries (&mut self) {
        if self.r[0] < -BOX_X  {
            self.v[0] = - self.elast * self.v[0];
            self.r[0] = -BOX_X * 0.99;
        }
        if self.r[0] > BOX_X  {
            self.v[0] = -self.elast * self.v[0];
            self.r[0] = BOX_X * 0.99;
        }
        if self.r[1] < -BOX_Y  {
            self.v[1] = -self.elast * self.v[1];
            self.r[1] = -BOX_Y * 0.99;
        }
        if self.r[1] > BOX_Y  {
            self.v[1] = -self.elast * self.v[1];
            self.r[1] = BOX_Y * 0.99;
        }
        if self.r[2] < -BOX_Z  {
            self.v[2] = -self.elast * self.v[2];
            self.r[2] = -BOX_Z * 0.99;
        }
        if self.r[2] > BOX_Z  {
            self.v[2] = -self.elast * self.v[2];
            self.r[2] = BOX_Z * 0.99;
        }
    }

    fn update(&mut self, positions: &Vec<[f64; 3]>, dt: f64) {
    // fn update(&mut self, a: &Vec<Particle>, dt: f64) {
        self.force_total();
        self.verlet(dt);
        self.check_boundaries();
    }

}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        let x = (args.width / 2) as f64;
        let y = (args.height / 2) as f64;
        let ref particles = self.particles;
        let ref bbox = self.bbox;

        self.gl.draw(
            args.viewport(), |c, gl| {
                clear(BACKGROUND, gl);

                let mut draw = |p: &Particle| {
                    let rad = p.projected_radius() * (args.height as f64)*0.5;
                    p.ellipse.draw(
                        [x+p.projected(p.r[0])-rad*0.5,
                         y-(p.projected(p.r[1])-rad*0.5),
                         rad, rad],
                        &c.draw_state, c.transform, gl);
                };

                for p in particles {
                    draw(p);
                }
                for p in bbox {
                    draw(p);
                }

        });
    }

    fn update(&mut self, args: &UpdateArgs) {
        let mut positions: Vec<[f64; 3]> = vec![];
        for p in self.particles.iter(){
            positions.push(p.r);
        }
        for p in self.particles.iter_mut(){
            p.update(&positions, args.dt);
        }
    }

}

fn main() {
    let opengl = OpenGL::V3_2;

    let window = Window::new(
        WindowSettings::new(
            "overflow", [1200, 800])
            .opengl(opengl)
            .exit_on_esc(true));

    let mut app = App {
        gl: GlGraphics::new(opengl),
        particles: vec![],
        bbox: vec![],
    };

    app.bbox.push(Particle::new(20.0, -BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  BOX_X, -BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0, -BOX_X,  BOX_Y,  BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0, -BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  BOX_X, -BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0, -BOX_X,  BOX_Y, -BOX_Z, 0.0, 0.0, 0.0));

    for i in  0..50 {
        app.particles.push(Particle::random_pos());
    }

    for e in window.events() {
        if let Some(r) = e.render_args() {app.render(&r);}
        if let Some(u) = e.update_args() {app.update(&u);}
    }
}
