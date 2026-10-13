use sims_field::Grid;
use sims_field::grid::{Grading, cell_of, centres, graded, nearest, pad, smooth, uniform};

#[test]
fn nearest_picks_the_closer_line_and_clamps_at_the_ends() {
    let line = [0.0, 1.0, 3.0, 7.0];
    assert_eq!(nearest(&line, -5.0), 0);
    assert_eq!(nearest(&line, 0.4), 0);
    assert_eq!(nearest(&line, 0.5), 0);
    assert_eq!(nearest(&line, 0.6), 1);
    assert_eq!(nearest(&line, 2.1), 2);
    assert_eq!(nearest(&line, 5.0), 2);
    assert_eq!(nearest(&line, 99.0), 3);
}

#[test]
fn cell_of_stays_inside_the_cells() {
    let line = [0.0, 1.0, 3.0, 7.0];
    assert_eq!(cell_of(&line, -1.0), 0);
    assert_eq!(cell_of(&line, 0.0), 0);
    assert_eq!(cell_of(&line, 1.0), 1);
    assert_eq!(cell_of(&line, 6.9), 2);
    assert_eq!(cell_of(&line, 7.0), 2);
    assert_eq!(cell_of(&line, 70.0), 2);
}

#[test]
fn a_uniform_grid_counts_its_cells() {
    let g = Grid::uniform([0.0; 3], [1.0, 2.0, 3.0], [2, 4, 6], 3);
    assert_eq!(g.dims(), [3, 5, 7]);
    assert_eq!(g.cells(), 48);
    assert_eq!(g.centres(0), vec![0.25, 0.75]);
    assert_eq!(uniform(0.0, 1.0, 4), vec![0.0, 0.25, 0.5, 0.75, 1.0]);
    assert_eq!(centres(&[0.0, 2.0, 6.0]), vec![1.0, 4.0]);
}

fn grade(fixed: &[f64], features: &[(f64, f64)], coarse: f64, ratio: f64, fill: f64) -> Vec<f64> {
    graded(fixed, &Grading { features, coarse: &|_| coarse, ratio, fill })
}

#[test]
fn graded_lines_keep_the_fixed_lines_and_grow_from_features() {
    let lines = grade(&[1.0, 0.0, 0.3, 0.3], &[(0.0, 0.01)], 0.1, 1.3, 0.0);
    for v in [0.0, 0.3, 1.0] {
        assert!(lines.iter().any(|x| (x - v).abs() < 1e-12), "{v} lost");
    }
    let d: Vec<f64> = lines.windows(2).map(|w| w[1] - w[0]).collect();
    assert!(d.iter().all(|d| *d > 0.0));
    assert!(d[0] < 0.012, "{}", d[0]);
    assert!(d.iter().all(|d| *d <= 0.1 + 1e-9));
    for w in d.windows(2) {
        assert!(w[1] / w[0] < 1.35, "{} then {}", w[0], w[1]);
    }
}

#[test]
fn the_fill_rule_drops_a_sliver_cell() {
    let plain = grade(&[0.0, 0.11], &[], 0.1, 1.3, 0.0);
    let filled = grade(&[0.0, 0.11], &[], 0.1, 1.3, 0.8);
    assert_eq!(plain.len(), 3);
    assert_eq!(filled, vec![0.0, 0.11]);
}

#[test]
fn smoothing_splits_a_jump() {
    let line = smooth(vec![0.0, 1.0, 2.0, 10.0], 1.5);
    for w in line.windows(3) {
        let (a, b) = (w[1] - w[0], w[2] - w[1]);
        assert!(b <= 1.5 * a + 1e-12 && a <= 1.5 * b + 1e-12, "{a} then {b}");
    }
    assert_eq!(*line.last().unwrap(), 10.0);
}

#[test]
fn padding_repeats_the_end_cells() {
    assert_eq!(pad(vec![0.0, 1.0, 3.0], 2), vec![-2.0, -1.0, 0.0, 1.0, 3.0, 5.0, 7.0]);
}
