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

struct Vec3d (f64, f64, f64);

pub struct Particle {
    r: Vec3d,
    v: Vec3d,
    a: Vec3d,
    a2: Vec3d,
    m: f64,
    rad: f64,
    ellipse: Ellipse,
}

// fn force_gravity (&mut p: Particle) {
//     p.
// }

// fn force_total (&mut p: Particle) {
//     p.
// }


impl Particle {
    fn new(rad: f64,
           rx: f64, ry: f64, rz: f64,
           vx: f64, vy: f64, vz: f64) -> Particle{
        Particle {
            r: Vec3d (rx, ry, rz),
            v: Vec3d (vx, vy, vz),
            a: Vec3d (0.0, 0.0, 0.0),
            a2: Vec3d (0.0, 0.0, 0.0),
            m: 1.0,
            rad: rad,
            ellipse: Ellipse::new([0.5, 0.5, 0.5, 1.0]),
        }
    }

    fn dist(&self) -> f64 {
        (self.r.0*self.r.0 + self.r.1*self.r.1 + self.r.2*self.r.2).sqrt()
    }

    fn projected_radius(&self) -> f64 {
        let fov = FOVY / 2.0 * PI / 180.0;
        let d = self.r.2 + SCREEN;
        1.0 / fov.tan() * self.rad / (d*d - self.rad*self.rad).sqrt()
    }

    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/(self.r.2 + SCREEN*2.0)
    }

    fn verlet(&mut self, dt: f64){
        self.r.0 = self.r.0 + self.v.0*dt + 0.5*self.a.0*dt*dt;
        self.r.1 = self.r.1 + self.v.1*dt + 0.5*self.a.1*dt*dt;
        self.r.2 = self.r.2 + self.v.2*dt + 0.5*self.a.2*dt*dt;

    }

    fn update(&mut self, dt: f64) {
        if self.r.0 < -BOX_X  {
            self.v.0 = -1.0 * self.v.0;
            self.r.0 = -BOX_X * 0.99;
        }
        if self.r.0 > BOX_X  {
            self.v.0 = -1.0 * self.v.0;
            self.r.0 = BOX_X * 0.99;
        }
        if self.r.1 < -BOX_Y  {
            self.v.1 = -1.0 * self.v.1;
            self.r.1 = -BOX_Y * 0.99;
        }
        if self.r.1 > BOX_Y  {
            self.v.1 = -1.0 * self.v.1;
            self.r.1 = BOX_Y * 0.99;
        }
        if self.r.2 < -BOX_Z  {
            self.v.2 = -1.0 * self.v.2;
            self.r.2 = -BOX_Z * 0.99;
        }
        if self.r.2 > BOX_Z  {
            self.v.2 = -1.0 * self.v.2;
            self.r.2 = BOX_Z * 0.99;
        }
        // self.verlet(dt);
    }

}


pub struct App {
    gl: GlGraphics,
    particles: Vec<Particle>,
    bbox: Vec<Particle>,
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
                    let xx = p.projected(p.r.0);
                    let yy = p.projected(p.r.1);
                    p.ellipse.draw(
                        [x+xx-rad*0.5, y-(yy-rad*0.5), rad, rad],
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
            self.particles[i].r.0 += self.particles[i].v.0 * args.dt;
            self.particles[i].r.1 += self.particles[i].v.1 * args.dt;
            self.particles[i].r.2 += self.particles[i].v.2 * args.dt;
            self.particles[i].update(args.dt);
            // let lum = ((PI/2.0 - (self.particles[i].dist()/200.).atan())
            //     /(PI/2.0)) as f32;
            // self.particles[i].ellipse.color = [lum, lum, lum, 1.0];
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

    app.particles.push(Particle::new(20.0, 0.0, 0.0, 0.0, 200.0, 200.0, 200.0));


    for e in window.events() {
        if let Some(r) = e.render_args() {app.render(&r);}
        if let Some(u) = e.update_args() {app.update(&u);}
    }
}
