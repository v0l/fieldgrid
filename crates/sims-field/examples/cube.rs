use sims_field::raster::{closed, solid_cells};
use sims_field::{Bvh, Grid};

fn main() {
    let p = |x: f64, y: f64, z: f64| [x, y, z];
    let quad = |a, b, c, d| [[a, b, c], [a, c, d]];
    let (lo, hi) = (0.2, 0.8);
    let mut cube = Vec::new();
    cube.extend(quad(p(lo, lo, lo), p(lo, hi, lo), p(hi, hi, lo), p(hi, lo, lo)));
    cube.extend(quad(p(lo, lo, hi), p(hi, lo, hi), p(hi, hi, hi), p(lo, hi, hi)));
    cube.extend(quad(p(lo, lo, lo), p(hi, lo, lo), p(hi, lo, hi), p(lo, lo, hi)));
    cube.extend(quad(p(lo, hi, lo), p(lo, hi, hi), p(hi, hi, hi), p(hi, hi, lo)));
    cube.extend(quad(p(lo, lo, lo), p(lo, lo, hi), p(lo, hi, hi), p(lo, hi, lo)));
    cube.extend(quad(p(hi, lo, lo), p(hi, hi, lo), p(hi, hi, hi), p(hi, lo, hi)));
    assert!(closed(&cube));

    let grid = Grid::uniform([0.0; 3], [1.0; 3], [20; 3], 0);
    let inside = solid_cells(&Bvh::new(cube), &grid);
    println!("{} of {} cells are inside", inside.iter().filter(|c| **c).count(), inside.len());
}
