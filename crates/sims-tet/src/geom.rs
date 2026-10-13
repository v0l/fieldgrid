pub type Point3 = cgmath::Point3<f64>;
pub type Vector3 = cgmath::Vector3<f64>;
pub type Vector2 = cgmath::Vector2<f64>;
pub type Matrix2 = cgmath::Matrix2<f64>;
pub type Matrix3 = cgmath::Matrix3<f64>;
pub use cgmath::prelude::*;

pub trait Surface {
    fn at(&self, uv: Vector2) -> Option<Point3>;
    fn nearest(&self, p: Point3, hint: Option<Vector2>) -> Option<Vector2>;
}

impl Surface for () {
    fn at(&self, _: Vector2) -> Option<Point3> {
        None
    }

    fn nearest(&self, _: Point3, _: Option<Vector2>) -> Option<Vector2> {
        None
    }
}
