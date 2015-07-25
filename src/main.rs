extern crate piston;
extern crate graphics;
extern crate glutin_window;
extern crate opengl_graphics;

use piston::window::WindowSettings;
use piston::event::*;
use glutin_window::GlutinWindow as Window;
use opengl_graphics::{ GlGraphics, OpenGL };
use graphics::*;
use std::f64::consts::PI;
use std::ops::{Sub, Add};


const BACKGROUND: [f32; 4] = [0.0, 0.0, 0.0, 1.0];
const FOVY: f64 = 60.0;
const BOX_X: f64 = 600.0;
const BOX_Y: f64 = 300.0;
const BOX_Z: f64 = 300.0;
const SCREEN: f64 = 600.0;

pub struct Particle {
    r: [f64; 3],      // location
    v: [f64; 3],      // velocity
    a: [f64; 3],      // acceleration
    ap: [f64; 3],     // acceleration from previous time-step
    m: f64,           // mass
    rad: f64,         // radius
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
            ellipse: Ellipse::new([0.5, 0.5, 0.5, 1.0]),
        }
    }

    fn dist(&self) -> f64 {
        (self.r[0]*self.r[0] + self.r[1]*self.r[1] + self.r[2]*self.r[2]).sqrt()
    }

    fn projected_radius(&self) -> f64 {
        let fov = FOVY / 2.0 * PI / 180.0;
        let d = self.r[2] + SCREEN;
        1.0 / fov.tan() * self.rad / (d*d - self.rad*self.rad).sqrt()
    }

    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/(self.r[2] + SCREEN*2.0)
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
            self.v[0] = -1.0 * self.v[0];
            self.r[0] = -BOX_X * 0.99;
        }
        if self.r[0] > BOX_X  {
            self.v[0] = -1.0 * self.v[0];
            self.r[0] = BOX_X * 0.99;
        }
        if self.r[1] < -BOX_Y  {
            self.v[1] = -1.0 * self.v[1];
            self.r[1] = -BOX_Y * 0.99;
        }
        if self.r[1] > BOX_Y  {
            self.v[1] = -1.0 * self.v[1];
            self.r[1] = BOX_Y * 0.99;
        }
        if self.r[2] < -BOX_Z  {
            self.v[2] = -1.0 * self.v[2];
            self.r[2] = -BOX_Z * 0.99;
        }
        if self.r[2] > BOX_Z  {
            self.v[2] = -1.0 * self.v[2];
            self.r[2] = BOX_Z * 0.99;
        }
    }

    fn update(&mut self, dt: f64) {
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
        for i in 0..self.particles.len(){
            self.particles[i].update(args.dt);
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

    app.particles.push(Particle::new(20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    for e in window.events() {
        if let Some(r) = e.render_args() {app.render(&r);}
        if let Some(u) = e.update_args() {app.update(&u);}
    }
}
