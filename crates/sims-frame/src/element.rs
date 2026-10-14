use crate::{Member, Model};

pub type M12 = [[f64; 12]; 12];
pub type V12 = [f64; 12];

pub struct Local {
    pub length: f64,
    pub rotation: [[f64; 3]; 3],
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = dot(a, a).sqrt();
    (n > 0.0).then(|| a.map(|v| v / n))
}

pub fn local(model: &Model, m: &Member) -> Result<Local, String> {
    let along = sub(model.nodes[m.b], model.nodes[m.a]);
    let length = dot(along, along).sqrt();
    let x = unit(along).ok_or("a member joins a node to itself")?;
    let up = if unit(cross(x, m.up)).is_some() {
        m.up
    } else if x[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let y = unit(sub(up, x.map(|v| v * dot(up, x)))).ok_or("a member's up runs along it")?;
    let z = cross(x, y);
    Ok(Local { length, rotation: [x, y, z] })
}

pub fn stiffness(m: &Member, l: f64) -> M12 {
    let (e, g) = (m.material.e, m.material.g);
    let s = &m.section;
    let mut k = [[0.0; 12]; 12];
    let mut put = |i: usize, j: usize, v: f64| {
        k[i][j] += v;
        if i != j {
            k[j][i] += v;
        }
    };
    let ea = e * s.area / l;
    put(0, 0, ea);
    put(0, 6, -ea);
    put(6, 6, ea);
    let gj = g * s.j / l;
    put(3, 3, gj);
    put(3, 9, -gj);
    put(9, 9, gj);
    let (z, l2, l3) = (e * s.iz, l * l, l * l * l);
    put(1, 1, 12.0 * z / l3);
    put(1, 5, 6.0 * z / l2);
    put(1, 7, -12.0 * z / l3);
    put(1, 11, 6.0 * z / l2);
    put(5, 5, 4.0 * z / l);
    put(5, 7, -6.0 * z / l2);
    put(5, 11, 2.0 * z / l);
    put(7, 7, 12.0 * z / l3);
    put(7, 11, -6.0 * z / l2);
    put(11, 11, 4.0 * z / l);
    let y = e * s.iy;
    put(2, 2, 12.0 * y / l3);
    put(2, 4, -6.0 * y / l2);
    put(2, 8, -12.0 * y / l3);
    put(2, 10, -6.0 * y / l2);
    put(4, 4, 4.0 * y / l);
    put(4, 8, 6.0 * y / l2);
    put(4, 10, 2.0 * y / l);
    put(8, 8, 12.0 * y / l3);
    put(8, 10, 6.0 * y / l2);
    put(10, 10, 4.0 * y / l);
    k
}

fn bending(target: &mut M12, dofs: [usize; 4], flip: f64, scale: f64, m: [[f64; 4]; 4]) {
    let sign = [1.0, flip, 1.0, flip];
    for i in 0..4 {
        for j in 0..4 {
            target[dofs[i]][dofs[j]] += scale * m[i][j] * sign[i] * sign[j];
        }
    }
}

pub fn mass(m: &Member, l: f64) -> M12 {
    let rho = m.material.density;
    let s = &m.section;
    let mut out = [[0.0; 12]; 12];
    let line = rho * s.area * l;
    let polar = rho * (s.iy + s.iz) * l;
    for (a, b, w) in [(0, 6, line), (3, 9, polar)] {
        out[a][a] += w / 3.0;
        out[b][b] += w / 3.0;
        out[a][b] += w / 6.0;
        out[b][a] += w / 6.0;
    }
    let l2 = l * l;
    let shape = [
        [156.0, 22.0 * l, 54.0, -13.0 * l],
        [22.0 * l, 4.0 * l2, 13.0 * l, -3.0 * l2],
        [54.0, 13.0 * l, 156.0, -22.0 * l],
        [-13.0 * l, -3.0 * l2, -22.0 * l, 4.0 * l2],
    ];
    bending(&mut out, [1, 5, 7, 11], 1.0, line / 420.0, shape);
    bending(&mut out, [2, 4, 8, 10], -1.0, line / 420.0, shape);
    out
}

pub fn geometric(axial: f64, l: f64) -> M12 {
    let mut out = [[0.0; 12]; 12];
    let l2 = l * l;
    let shape = [
        [36.0, 3.0 * l, -36.0, 3.0 * l],
        [3.0 * l, 4.0 * l2, -3.0 * l, -l2],
        [-36.0, -3.0 * l, 36.0, -3.0 * l],
        [3.0 * l, -l2, -3.0 * l, 4.0 * l2],
    ];
    bending(&mut out, [1, 5, 7, 11], 1.0, axial / (30.0 * l), shape);
    bending(&mut out, [2, 4, 8, 10], -1.0, axial / (30.0 * l), shape);
    out
}

pub fn uniform(w: [f64; 3], l: f64) -> V12 {
    let mut q = [0.0; 12];
    q[0] = w[0] * l / 2.0;
    q[6] = w[0] * l / 2.0;
    q[1] = w[1] * l / 2.0;
    q[7] = w[1] * l / 2.0;
    q[5] = w[1] * l * l / 12.0;
    q[11] = -w[1] * l * l / 12.0;
    q[2] = w[2] * l / 2.0;
    q[8] = w[2] * l / 2.0;
    q[4] = -w[2] * l * l / 12.0;
    q[10] = w[2] * l * l / 12.0;
    q
}

pub fn released(m: &Member) -> Vec<usize> {
    let mut out = Vec::new();
    if m.pinned[0] {
        out.extend([4, 5]);
    }
    if m.pinned[1] {
        out.extend([10, 11]);
    }
    out
}

pub struct Condensed {
    kept: Vec<usize>,
    map: Vec<V12>,
}

fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut row = r.clone();
            row.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 }));
            row
        })
        .collect();
    for c in 0..n {
        let p = (c..n).max_by(|x, y| m[*x][c].abs().total_cmp(&m[*y][c].abs()))?;
        if m[p][c].abs() < 1e-300 {
            return None;
        }
        m.swap(c, p);
        let d = m[c][c];
        m[c].iter_mut().for_each(|v| *v /= d);
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                let pivot = m[c].clone();
                m[r].iter_mut().zip(&pivot).for_each(|(v, p)| *v -= f * p);
            }
        }
    }
    Some(m.into_iter().map(|r| r[n..].to_vec()).collect())
}

impl Condensed {
    pub fn new(k: &M12, released: &[usize]) -> Result<Condensed, String> {
        let kept: Vec<usize> = (0..12).filter(|i| !released.contains(i)).collect();
        let mut map = [[0.0; 12]; 12];
        for &i in &kept {
            map[i][i] = 1.0;
        }
        if !released.is_empty() {
            let krr: Vec<Vec<f64>> =
                released.iter().map(|&i| released.iter().map(|&j| k[i][j]).collect()).collect();
            let inv = invert(&krr).ok_or("a member is released at both ends about one axis")?;
            for (a, &r) in released.iter().enumerate() {
                for &c in &kept {
                    map[r][c] =
                        -(0..released.len()).map(|b| inv[a][b] * k[released[b]][c]).sum::<f64>();
                }
            }
        }
        Ok(Condensed { kept, map: map.to_vec() })
    }

    pub fn matrix(&self, a: &M12) -> M12 {
        let mut out = [[0.0; 12]; 12];
        for &i in &self.kept {
            for &j in &self.kept {
                let mut s = 0.0;
                for p in 0..12 {
                    if self.map[p][i] == 0.0 {
                        continue;
                    }
                    for q in 0..12 {
                        s += self.map[p][i] * a[p][q] * self.map[q][j];
                    }
                }
                out[i][j] = s;
            }
        }
        out
    }

    pub fn vector(&self, v: &V12) -> V12 {
        let mut out = [0.0; 12];
        for &i in &self.kept {
            out[i] = (0..12).map(|p| self.map[p][i] * v[p]).sum();
        }
        out
    }
}

pub fn transform(r: &[[f64; 3]; 3]) -> M12 {
    let mut t = [[0.0; 12]; 12];
    for b in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                t[3 * b + i][3 * b + j] = r[i][j];
            }
        }
    }
    t
}

pub fn to_global(t: &M12, a: &M12) -> M12 {
    let mut ta = [[0.0; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            ta[i][j] = (0..12).map(|p| t[p][i] * a[p][j]).sum();
        }
    }
    let mut out = [[0.0; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            out[i][j] = (0..12).map(|p| ta[i][p] * t[p][j]).sum();
        }
    }
    out
}

pub fn vector_to_global(t: &M12, v: &V12) -> V12 {
    std::array::from_fn(|i| (0..12).map(|p| t[p][i] * v[p]).sum())
}

pub fn vector_to_local(t: &M12, v: &V12) -> V12 {
    std::array::from_fn(|i| (0..12).map(|p| t[i][p] * v[p]).sum())
}

pub fn rotate(r: &[[f64; 3]; 3], w: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| dot(r[i], w))
}
