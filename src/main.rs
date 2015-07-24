extern crate nalgebra as na;
use na::{Vec3, Rot3, Rotation};


fn main() {

    let     a = Vec3::new(1.0f64, 1.0, 1.0);
    let mut b = Rot3::new(na::zero());

    b.append_rotation_mut(&a);

    assert!(na::approx_eq(&na::rotation(&b), &a));
    println!("{}", a.x);

}
