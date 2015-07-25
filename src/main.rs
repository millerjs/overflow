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

struct Vec3d {
    x: f64,
    y: f64,
    z: f64,
}

pub struct Particle {
    r: Vec3d,
    v: Vec3d,
    a: Vec3d,
    a2: Vec3d,
    m: f64,
    rad: f64,
    ellipse: Ellipse,
}

impl Add for Vec3d {
    type Output = Vec3d;
    fn add(self, other: Vec3d) -> Vec3d {
        Vec3d {x: self.x + other.x, y: self.y + other.y, z: self.z + other.z}
    }
}

impl Sub for Vec3d {
    type Output = Vec3d;
    fn sub(self, other: Vec3d) -> Vec3d {
        Vec3d {x: self.x - other.x, y: self.y - other.y, z: self.z - other.z}
    }
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
            r: Vec3d {x: rx, y: ry, z: rz},
            v: Vec3d {x: vx, y: vy, z: vz},
            a: Vec3d {x: 0.0, y: 0.0, z: 0.0},
            a2: Vec3d {x: 0.0, y: 0.0, z: 0.0},
            m: 1.0,
            rad: rad,
            ellipse: Ellipse::new([0.5, 0.5, 0.5, 1.0]),
        }
    }

    fn dist(&self) -> f64 {
        (self.r.x*self.r.x + self.r.y*self.r.y + self.r.z*self.r.z).sqrt()
    }

    fn projected_radius(&self) -> f64 {
        let fov = FOVY / 2.0 * PI / 180.0;
        let d = self.r.z + SCREEN;
        1.0 / fov.tan() * self.rad / (d*d - self.rad*self.rad).sqrt()
    }

    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/(self.r.z + SCREEN*2.0)
    }

    fn verlet(&mut self, dt: f64){
        self.r.x = self.r.x + self.v.x*dt + 0.5*self.a.x*dt*dt;
        self.r.y = self.r.y + self.v.y*dt + 0.5*self.a.y*dt*dt;
        self.r.z = self.r.z + self.v.z*dt + 0.5*self.a.z*dt*dt;

    }

    fn update(&mut self, dt: f64) {
        if self.r.x < -BOX_X  {
            self.v.x = -1.0 * self.v.x;
            self.r.x = -BOX_X * 0.99;
        }
        if self.r.x > BOX_X  {
            self.v.x = -1.0 * self.v.x;
            self.r.x = BOX_X * 0.99;
        }
        if self.r.y < -BOX_Y  {
            self.v.y = -1.0 * self.v.y;
            self.r.y = -BOX_Y * 0.99;
        }
        if self.r.y > BOX_Y  {
            self.v.y = -1.0 * self.v.y;
            self.r.y = BOX_Y * 0.99;
        }
        if self.r.z < -BOX_Z  {
            self.v.z = -1.0 * self.v.z;
            self.r.z = -BOX_Z * 0.99;
        }
        if self.r.z > BOX_Z  {
            self.v.z = -1.0 * self.v.z;
            self.r.z = BOX_Z * 0.99;
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
                    let xx = p.projected(p.r.x);
                    let yy = p.projected(p.r.y);
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
            self.particles[i].r.x += self.particles[i].v.x * args.dt;
            self.particles[i].r.y += self.particles[i].v.y * args.dt;
            self.particles[i].r.z += self.particles[i].v.z * args.dt;
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
