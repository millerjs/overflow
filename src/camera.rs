extern crate rand;
extern crate image;

// Usage statements
use std::f64::consts::PI;
use std::thread;
use image::{ImageBuffer, Rgb};
use std::cmp::{min, max};

use linalg::Matrix;
use linalg::Vec3d;
use overflow::Camera;
use overflow::Particle;


#[allow(dead_code)]
impl Camera {

    pub fn new(width: i32, height: i32, r: Vec3d) -> Camera {
        println!("Creating image buffer of length {} KB", width*height/1000);

        Camera {
            pixels: vec![[0; 4]; width as usize * height as usize],
            render_step: 0,
            width: width,
            height: height,
            r: r,
            screen: 3.0,
            fovd: 60.0,
            theta: Vec3d::new(0.0, 0.0, 0.0),
        }
    }

    /// Scale the apparent radius of the ellipse for perspective
    pub fn projected_radius(&self, r: Vec3d, rad: f64) -> f64 {
        let fov = self.fovd / 2.0 * PI / 180.0;
        let d = (r - self.r).norm();
        (2e-1 / fov.tan() * rad /
            (d*d - rad*rad).sqrt() * (self.width) as f64).abs()
    }

    /// Get the x, y location as projected on the screen
    pub fn projected(&self, r: Vec3d) -> [i32; 2] {
        let d = &Matrix::projection(self.theta) * (r-self.r);
        let bx = (d[0] * self.width  as f64)/(d[2]*self.screen);
        let by = (d[1] * self.width as f64)/(d[2]*self.screen);
        [self.width/2 + bx as i32, self.height/2 - by as i32]
    }

    pub fn set_pixel(&mut self, x: i32, y: i32, pixel: [f32; 4]) {
        if x < self.width && x > 0 && y < self.height && y > 0 {
            for i in 0..4 {
                let idx = (self.height * x + y) as usize;
                self.pixels[idx][i] = (pixel[i]*255.0) as u8;
            }
        }
    }

    pub fn get_pixel(&mut self, x: i32, y: i32) -> [u8; 4]{
        if x < self.width && x > 0 && y < self.height && y > 0 {
            let idx = (self.height * x + y) as usize;
            self.pixels[idx]
        } else {
            [0; 4]
        }
    }

    pub fn is_visible(&self, r: Vec3d) -> bool {
        let d = &Matrix::projection(self.theta) * (r-self.r);
        d[2] > 0.0
    }

    pub fn draw_particle(&mut self, p: &Particle, pixel: [f32; 4]) {
        if !self.is_visible(p.r){ return };

        let translated = self.projected(p.r);
        let (x0, y0) = (translated[0], translated[1]);
        let mut x = self.projected_radius(p.r, p.rad) as i32;
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

    pub fn draw_line(&mut self, _p0: Vec3d, _p1: Vec3d, pixel: [f32; 4]) {
        if !self.is_visible(_p0){ return };
        if !self.is_visible(_p1){ return };

        let p0 = self.projected(_p0);
        let p1 = self.projected(_p1);
        let (x0, y0, x1, y1) = (p0[0], p0[1], p1[0], p1[1]);
        let (dx, dy) = (x1 - x0, y1 - y0);
        if dx < dy {
            for x in (min(x0, x1)..max(x0, x1)) {
                self.set_pixel(x, y0+(x-x0)*dy/dx, pixel);
            }
        } else {
            for y in (min(y0, y1)..max(y0, y1)) {
                self.set_pixel(x0+(y-y0)*dx/dy, y, pixel);
            }
        }
    }


    pub fn clear(&mut self, color: [f32; 4]) {
        for i in 0..(self.width*self.height) as usize {
            for j in 0..4 {
                self.pixels[i][j] = (color[j]*255 as f32) as u8;
            }
        }
    }

    pub fn save(&mut self) -> thread::JoinHandle<()> {
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

        self.render_step += 1;

        let postfix = self.render_step;
        thread::spawn(move || {
            let path = format!("output/test_{:08}.png", postfix);
            imgbuf.save(&*path).unwrap();
        })
    }

    pub fn rotate(&mut self, x: f64, y: f64, z: f64) {
        self.theta[0] += x;
        self.theta[1] += y;
        self.theta[2] += z;
    }

    pub fn revolve_y(&mut self, theta: f64) {
        self.theta[1] -= theta;
        self.r = &Matrix::rot_y(-theta)*self.r;
    }

}
