use std::ops;

#[derive(Debug)]
pub struct Matrix<T> {
    pub n: usize,
    pub m: usize,
    pub data: Vec<T>,
}

#[derive(Debug, Copy, Clone)]
pub struct Vec3d {
    pub data: [f64; 3],
}

// #########################################################################
// Matrices
// #########################################################################

#[allow(dead_code)]
impl Matrix<f64> {

    pub fn new(n: usize, m: usize, val: f64) -> Matrix<f64> {
        Matrix {n: n, m: m, data: vec![val; n*m]}
    }

    pub fn rot_x(theta: f64) -> Matrix<f64> {
        Matrix::from_data(3, 3, vec![
            1.0,     0.0,          0.0,
            0.0, theta.cos(), -theta.sin(),
            0.0, theta.sin(),  theta.cos()])
    }

    pub fn rot_y(theta: f64) -> Matrix<f64> {
        Matrix::from_data(3, 3, vec![
             theta.cos(), 0.0, theta.sin(),
                0.0,      1.0,     0.0,
            -theta.sin(), 0.0, theta.cos()])
    }

    pub fn rot_z(theta: f64) -> Matrix<f64> {
        Matrix::from_data(3, 3, vec![
            theta.cos(), -theta.sin(), 0.0,
            theta.sin(),  theta.cos(), 0.0,
                0.0,          0.0,     1.0])
    }

    pub fn from(array: Vec<Vec<f64>>) -> Matrix<f64> {
        let m = array.len();
        let n = array[0].len();
        let mut new = Matrix::zero(n, m);
        for i in 0..n {
            for j in 0..m {
                new[(i, j)] = array[j][i];
            }
        }
        new
    }

    pub fn from_data(n: usize, m: usize, array: Vec<f64>) -> Matrix<f64> {
        let mut new = Matrix::zero(n, m);
        for i in 0..n*m { new.data[i] = array[i] };
        new
    }

    pub fn zero(n: usize, m: usize) -> Matrix<f64> {
        Matrix::new(n, m, 0.0)
    }

    pub fn identity(n: usize) -> Matrix<f64> {
        let mut mat = Matrix::new(n, n, 0.0);
        for i in 0..n { mat[(i, i)] = 1.0 }
        mat
    }

    pub fn t(&self) -> Matrix<f64> {
        let mut mat = Matrix::zero(self.m, self.n);
        for i in 0..self.m {
            for j in 0..self.n {
                mat[(i, j)] = self[(j, i)];
            }
        }
        mat
    }

    pub fn to_vec3d(&self) -> Vec3d {
        assert_eq!(self.n, 1);
        assert_eq!(self.m, 3);
        Vec3d::new(self.data[0], self.data[1], self.data[2])
    }

    pub fn copy(&self) -> Matrix<f64> {
        let mut mat = Matrix::zero(self.n, self.m);
        for i in 0..self.n {
            for j in 0..self.m {
                mat[(i, j)] = self[(i, j)];
            }
        }
        mat
    }

    pub fn scale(&mut self, s: f64) {
        for i in 0..self.n*self.m {
            self.data[i] *= s;
        };
    }

    pub fn add(&mut self, s: f64) {
        for i in 0..self.n*self.m {
            self.data[i] += s;
        };
    }

    pub fn sum(&mut self) -> f64 {
        (0..self.n*self.m).fold(0.0, |s, i| s + self.data[i])
    }

    pub fn print(&self) {
        print!("[");
        for j in 0..self.m {
            print!("{}", if j==0 {"["} else {" ["});
            for i in 0..self.n {
                if i < self.n-1 {
                    print!("{}, ", self[(i, j)]);
                } else {
                    print!("{}", self[(i, j)]);
                }
            }
            print!("{}", if j<self.m-1 {"]\n"} else {"]"});
        }
        println!("]");
    }
}

impl ops::Index<(usize, usize)> for Matrix<f64> {
    type Output = f64;
    fn index<'a>(&'a self, ij: (usize, usize)) -> &'a f64 {
        &self.data[ij.1*self.n + ij.0]
    }
}

impl ops::IndexMut<(usize, usize)> for Matrix<f64> {
    fn index_mut<'a>(&'a mut self, ij: (usize, usize)) -> &'a mut f64 {
        &mut self.data[ij.1*self.n + ij.0]
    }
}

impl ops::Add for Matrix<f64> {
    type Output = Matrix<f64>;
    fn add(self, other: Matrix<f64>) -> Matrix<f64> {
        assert_eq!(self.n, other.n);
        assert_eq!(self.m, other.m);
        let mut new = Matrix::zero(self.n, self.m);
        for i in 0..self.n {
            for j in 0..self.m {
                new[(i, j)] = self[(i, j)] + other[(i, j)]
            }
        }
        new
    }
}

impl<'a> ops::Mul<&'a Matrix<f64>> for &'a Matrix<f64> {
    type Output = Matrix<f64>;
    fn mul(self, other: &Matrix<f64>) -> Matrix<f64> {
        assert_eq!(self.n, other.m);
        let mut new = Matrix::zero(other.n, self.m);
        for i in 0..other.n {
            for j in 0..self.m {
                for k in 0..other.m {
                    new[(i, j)] += other[(i, k)] * self[(k, j)]
                }
            }
        }
        new
    }
}

impl<'a> ops::Mul<Vec3d> for &'a Matrix<f64> {
    type Output = Vec3d;
    fn mul(self, v: Vec3d) -> Vec3d {
        let other = v.to_matrix();
        let mut new = Matrix::zero(other.n, self.m);
        for i in 0..other.n {
            for j in 0..self.m {
                for k in 0..other.m {
                    new[(i, j)] += other[(i, k)] * self[(k, j)]
                }
            }
        }
        new.to_vec3d()
    }
}

impl<'a> ops::Mul<f64> for &'a Matrix<f64> {
    type Output = Matrix<f64>;
    fn mul(self, scale: f64) -> Matrix<f64> {
        let mut new = self.copy();
        new.scale(scale);
        new
    }
}

impl<'a> ops::Mul<&'a Matrix<f64>> for f64 {
    type Output = Matrix<f64>;
    fn mul(self, m: &Matrix<f64>) -> Matrix<f64> {
        let mut new = m.copy();
        new.scale(self);
        new
    }
}


// #########################################################################
// Vectors
// #########################################################################

#[allow(dead_code)]
impl Vec3d {

    pub fn new(x: f64, y: f64, z: f64) -> Vec3d {
        Vec3d { data: [x, y, z]}
    }

    pub fn from(vals: [f64; 3]) -> Vec3d {
        Vec3d {data: vals}
    }

    pub fn zero() -> Vec3d {
        Vec3d {data: [0.0; 3]}
    }

    pub fn copy(&self) -> Vec3d {
        Vec3d {data: self.data}
    }

    pub fn scale(&mut self, s: f64) {
        for i in 0..3 { self.data[i] *= s };
    }

    pub fn add(&mut self, other: Vec3d) {
        for i in 0..3 { self.data[i] += other.data[i] };
    }

    pub fn sub(&mut self, other: Vec3d) {
        for i in 0..3 { self.data[i] -= other.data[i] };
    }

    pub fn sum(&self) -> f64 {
        self.data[0] + self.data[1] + self.data[2]
    }

    pub fn to_matrix(&self) -> Matrix<f64> {
        Matrix::from_data(1, 3, self.data.to_vec())
    }

    pub fn norm(&self) -> f64 {
        (self[0]*self[0] + self[1]*self[1] + self[2]*self[2]).sqrt()
    }

    pub fn unit(&self) -> Vec3d {
        let norm = self.norm();
        Vec3d::new(self[0]/norm, self[1]/norm, self[2]/norm)
    }

    // pub fn rot_x(&self, theta: f64) -> Vec3d {
    //     let r = Matrix::rot_x(theta);
    // }

}

impl ops::Mul<f64> for Vec3d {
    type Output = Vec3d;
    fn mul(self, scale: f64) -> Vec3d {
        let mut new = self.copy();
        new.scale(scale);
        new
    }
}

impl ops::BitXor<Vec3d> for Vec3d {
    type Output = Vec3d;
    fn bitxor(self, other: Vec3d) -> Vec3d {
        Vec3d::new(self[1]*other[2] - self[2]*other[1],
                   self[2]*other[0] - self[0]*other[2],
                   self[0]*other[1] - self[1]*other[0])
    }
}

impl ops::Mul<Vec3d> for f64 {
    type Output = Vec3d;
    fn mul(self, v: Vec3d) -> Vec3d {
        let mut new = v.copy();
        new.scale(self);
        new
    }
}

impl ops::Add<Vec3d> for Vec3d {
    type Output = Vec3d;
    fn add(self, other: Vec3d) -> Vec3d {
        Vec3d::new(self[0]+other[0], self[1]+other[1], self[2]+other[2])
    }
}

impl ops::Sub<Vec3d> for Vec3d {
    type Output = Vec3d;
    fn sub(self, other: Vec3d) -> Vec3d {
        Vec3d::new(self[0]-other[0], self[1]-other[1], self[2]-other[2])
    }
}

impl ops::Mul<Vec3d> for Vec3d {
    type Output = f64;
    fn mul(self, other: Vec3d) -> f64 {
        self[0]*other[0] + self[1]*other[1] + self[2]*other[2]
    }
}


impl ops::Index<usize> for Vec3d {
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64 {
        &self.data[i]
    }
}

impl ops::IndexMut<usize> for Vec3d {
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut f64 {
        &mut self.data[i]
    }
}
