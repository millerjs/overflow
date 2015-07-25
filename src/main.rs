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

const BACKGROUND: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

struct Vec3d {
    x: f64,
    y: f64,
    z: f64,
}

pub struct Particle {
    r: Vec3d,
    v: Vec3d
}

impl Particle {
    fn new(rx: f64, ry: f64, rz: f64, vx: f64, vy: f64, vz: f64) -> Particle{
        Particle {
            r: Vec3d {x: rx, y: ry, z: rz},
            v: Vec3d {x: vx, y: vy, z: vz}
        }
    }
}


pub struct App {
    gl: GlGraphics,
    particles: Vec<Particle>,
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        let x = (args.width / 3) as f64;
        let y = (args.height / 2) as f64;

        self.gl.draw(
            args.viewport(), |c, gl| {
                clear(BACKGROUND, gl);

                let transform = c.transform.trans(100.0, 100.0);

                let _ellipse = Ellipse::new([0.5, 0.5, 0.5, 1.0]);
                _ellipse.draw([x, y, 20.0, 20.0], &c.draw_state,
                              c.transform, gl);


        });
    }

    fn update(&mut self, args: &UpdateArgs) {
        // self.square.x += 10.0 * args.dt;
        // self.square.y += 10.0 * args.dt;
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

    app.particles.push(Particle::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

    for e in window.events() {
        if let Some(r) = e.render_args() {app.render(&r);}
        if let Some(u) = e.update_args() {app.update(&u);}
    }
}
