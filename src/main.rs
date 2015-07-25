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

const BLACK: [f32; 4] = [0.0, 0.0, 0.0, 1.0];
const RED:   [f32; 4] = [0.3, 0.3, 0.3, 1.0];

pub struct Square {
    x: f64,
    y: f64,
    rot: f64,
}


pub struct App {
    gl: GlGraphics,
    square: Square,
}

impl App {
    fn render(&mut self, args: &RenderArgs) {
        use graphics::*;

        let square = rectangle::square(0.0, 0.0, 50.0);
        let x = (args.width / 3) as f64 + self.square.x;
        let y = (args.height / 2) as f64 + self.square.y;

        self.gl.draw(
            args.viewport(), |c, gl| {
                clear(BLACK, gl);
                let transform = c.transform.trans(x, y);
                rectangle(RED, square, transform, gl);
        });
    }

    fn update(&mut self, args: &UpdateArgs) {
        self.square.x += 10.0 * args.dt;
        self.square.y += 10.0 * args.dt;
    }

}

fn main() {
    let opengl = OpenGL::V3_2;

    // Create an Glutin window.
    let window = Window::new(
        WindowSettings::new(
            "spinning-square", [200, 200])
            .opengl(opengl)
            .exit_on_esc(true));

    // Create a new game and run it.
    let mut app = App {
        gl: GlGraphics::new(opengl),
        square: Square{x: 0.0, y: 0.0, rot: 0.0}
    };

    for e in window.events() {
        if let Some(r) = e.render_args() {
            app.render(&r);
        }

        if let Some(u) = e.update_args() {
            app.update(&u);
        }
    }
}
