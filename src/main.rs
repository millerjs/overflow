extern crate nalgebra as na;

// use na::{Vec3, Rot3, Rotation};

// fn main() {
//     let     a = Vec3::new(1.0f64, 1.0, 1.0);
//     let mut b = Rot3::new(na::zero());

//     b.append_rotation_mut(&a);

//     assert!(na::approx_eq(&na::rotation(&b), &a));
//     println!("{}", a.x);

// }

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

const BACKGROUND: [f32; 4] = [0.0, 0.0, 0.0, 1.0];
const D2PIXEL: f64 = 100.0;

struct Vec3d {
    x: f64,
    y: f64,
    z: f64,
}

pub struct Particle {
    r: Vec3d,
    v: Vec3d,
    rad: f64,
    ellipse: Ellipse,
}

impl Particle {
    fn new(rad: f64,
           rx: f64, ry: f64, rz: f64,
           vx: f64, vy: f64, vz: f64) -> Particle{
        Particle {
            r: Vec3d {x: rx, y: ry, z: rz},
            v: Vec3d {x: vx, y: vy, z: vz},
            rad: rad,
            ellipse: Ellipse::new([0.5, 0.5, 0.5, 1.0]),
        }
    }

    fn dist(&self) -> f64 {
        (self.r.x*self.r.x + self.r.y*self.r.y + self.r.z*self.r.z).sqrt()
    }

    fn projectedRadius(&self) -> f64 {
        let fovy = 60.0;
        let fov = fovy / 2.0 * PI / 180.0;
        let d = self.dist();
        1.0 / fov.tan() * self.rad / (d*d + self.rad*self.rad).sqrt()
    }

}


pub struct App {
    gl: GlGraphics,
    particles: Vec<Particle>,
}

impl App {
    fn render(&mut self, args: &RenderArgs) {

        let x = (args.width / 3) as f64;
        let y = (args.height / 2) as f64;
        let ref particles = self.particles;

        self.gl.draw(
            args.viewport(), |c, gl| {
                clear(BACKGROUND, gl);

                for p in particles.iter() {
                    let rad = p.projectedRadius();
                    // println!("radius: {}", rad);
                    p.ellipse.draw(
                        [p.r.x, p.r.y, rad, rad],
                        &c.draw_state, c.transform, gl);
                }

        });
    }

    fn update(&mut, self, args: &UpdateArgs) {
        for p in self.particles.iter() {
            p.r.x += 10.0 * args.dt;
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
        particles: vec![]
    };

    app.particles.push(Particle::new(20.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0));
    app.particles.push(Particle::new(20.0, 100.0, 200.0, 0.0, 0.0, 0.0, 0.0));

    for e in window.events() {
        if let Some(r) = e.render_args() {app.render(&r);}
        if let Some(u) = e.update_args() {app.update(&u);}
    }
}
