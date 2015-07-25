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
const SCREEN: f64 = 20.0;

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

    fn projected_radius(&self) -> f64 {
        let fovy = 60.0;
        let fov = fovy / 2.0 * PI / 180.0;
        let d = self.dist();
        1.0 / fov.tan() * self.rad / (d*d - self.rad*self.rad).sqrt()
    }

    fn projected(&self, r: f64) -> f64 {
        SCREEN*r/self.r.z
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
                        [x+xx-rad*0.5, y+yy-rad*0.5, rad, rad],
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
            let lum = ((PI/2.0 - (self.particles[i].dist()/200.).atan())
                /(PI/2.0)) as f32;
            self.particles[i].ellipse.color = [lum, lum, lum, 1.0];
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


    app.bbox.push(Particle::new(20.0, -200.0, -200.0, 200.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  200.0, -200.0, 200.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  200.0,  200.0, 200.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0, -200.0,  200.0, 200.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0, -200.0, -200.0, 20.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  200.0, -200.0, 20.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0,  200.0,  200.0, 20.0, 0.0, 0.0, 0.0));
    app.bbox.push(Particle::new(20.0, -200.0,  200.0, 20.0, 0.0, 0.0, 0.0));

    // app.particles.push(Particle::new(20.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0));


    for e in window.events() {
        if let Some(r) = e.render_args() {app.render(&r);}
        if let Some(u) = e.update_args() {app.update(&u);}
    }
}
