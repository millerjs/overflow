extern crate image;
extern crate rand;
extern crate threadpool;

mod linalg;
mod overflow;
mod camera;
mod domain;
mod particle;

use std::thread;
use linalg::Vec3d;
use overflow::Camera;
use overflow::Domain;


fn main() {

    let max_time = 10.0;
    let mut step = 0;

    let mut domain = Domain {
        camera: Camera::new(300, 200, Vec3d::new(0.0, 0.0, -1000.0)),
        particles: vec![],
        bbox: vec![],
        dt: 1.0e-2,
        t: 0.0,
    };

    println!("Created image buffer of length {} KB",
             domain.camera.width*domain.camera.height/1000);

    domain.add_particles_random(25.0, 1000);

    while domain.t < max_time {
        domain.update();
        if step % 5 == 0 {
            domain.print_state();
            domain.render();
            domain.camera.revolve_y(0.05);
        }
        step += 1;
    }

    thread::sleep_ms(2000);

}
