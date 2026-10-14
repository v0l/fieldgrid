mod air;
mod decay;

pub use air::Air;
pub use decay::{MODULATION, Response};

use rayon::prelude::*;
use sims_field::Bvh;

pub const BANDS: [f64; 7] = [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0];

pub type Triangle = [[f64; 3]; 3];

#[derive(Clone, Debug, PartialEq)]
pub struct Surface {
    pub triangles: Vec<Triangle>,
    pub absorption: [f64; BANDS.len()],
    pub scattering: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Room {
    pub surfaces: Vec<Surface>,
    pub air: Air,
    pub air_absorption: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rays {
    pub count: usize,
    pub seconds: f64,
    pub bin: f64,
    pub radius: f64,
    pub seed: u64,
}

impl Default for Rays {
    fn default() -> Rays {
        Rays { count: 20_000, seconds: 2.0, bin: 0.001, radius: 0.5, seed: 1 }
    }
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

fn scale(a: [f64; 3], k: f64) -> [f64; 3] {
    a.map(|v| v * k)
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn normal(t: &Triangle) -> [f64; 3] {
    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    scale(n, 1.0 / dot(n, n).sqrt())
}

fn area(t: &Triangle) -> f64 {
    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    0.5 * dot(n, n).sqrt()
}

impl Room {
    pub fn volume(&self) -> f64 {
        self.surfaces
            .iter()
            .flat_map(|s| &s.triangles)
            .map(|t| dot(t[0], cross(t[1], t[2])) / 6.0)
            .sum::<f64>()
            .abs()
    }

    pub fn area(&self) -> f64 {
        self.surfaces.iter().flat_map(|s| &s.triangles).map(area).sum()
    }

    pub fn absorption_area(&self, band: usize) -> f64 {
        self.surfaces
            .iter()
            .map(|s| s.absorption[band] * s.triangles.iter().map(area).sum::<f64>())
            .sum()
    }

    fn air_term(&self, band: usize) -> f64 {
        if self.air_absorption {
            4.0 * self.air.nepers_per_metre(BANDS[band]) * self.volume()
        } else {
            0.0
        }
    }

    pub fn sabine(&self, band: usize) -> f64 {
        24.0 * std::f64::consts::LN_10 / self.air.speed() * self.volume()
            / (self.absorption_area(band) + self.air_term(band))
    }

    pub fn eyring(&self, band: usize) -> f64 {
        let s = self.area();
        let mean = self.absorption_area(band) / s;
        24.0 * std::f64::consts::LN_10 / self.air.speed() * self.volume()
            / (-s * (1.0 - mean).ln() + self.air_term(band))
    }
}

struct Random(u64);

impl Random {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }

    fn sphere(&mut self) -> [f64; 3] {
        let z = 2.0 * self.next() - 1.0;
        let phi = 2.0 * std::f64::consts::PI * self.next();
        let r = (1.0 - z * z).sqrt();
        [r * phi.cos(), r * phi.sin(), z]
    }

    fn lambert(&mut self, n: [f64; 3]) -> [f64; 3] {
        let (u1, u2) = (self.next(), self.next());
        let (r, phi) = (u1.sqrt(), 2.0 * std::f64::consts::PI * u2);
        let helper = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let a = cross(n, helper);
        let a = scale(a, 1.0 / dot(a, a).sqrt());
        let b = cross(n, a);
        add(add(scale(a, r * phi.cos()), scale(b, r * phi.sin())), scale(n, (1.0 - u1).sqrt()))
    }
}

fn chord(
    origin: [f64; 3],
    direction: [f64; 3],
    length: f64,
    centre: [f64; 3],
    radius: f64,
) -> Option<(f64, f64)> {
    let oc = sub(origin, centre);
    let b = dot(oc, direction);
    let c = dot(oc, oc) - radius * radius;
    let disc = b * b - c;
    if disc <= 0.0 {
        return None;
    }
    let root = disc.sqrt();
    let (t0, t1) = ((-b - root).max(0.0), (-b + root).min(length));
    (t1 > t0).then_some((t0, t1 - t0))
}

pub fn trace(
    room: &Room,
    source: [f64; 3],
    receivers: &[[f64; 3]],
    rays: &Rays,
) -> Result<Vec<Response>, String> {
    let stage = sims_progress::stage("tracing rays", Some(rays.count), "rays");
    let mut triangles = Vec::new();
    let mut owner = Vec::new();
    for (k, s) in room.surfaces.iter().enumerate() {
        triangles.extend(s.triangles.iter().copied());
        owner.extend(std::iter::repeat_n(k, s.triangles.len()));
    }
    if triangles.is_empty() {
        return Err("the room has no surfaces".into());
    }
    let normals: Vec<[f64; 3]> = triangles.iter().map(normal).collect();
    let bvh = Bvh::new(triangles);
    let c = room.air.speed();
    let bins = (rays.seconds / rays.bin).ceil() as usize;
    let air: [f64; BANDS.len()] = std::array::from_fn(|b| {
        if room.air_absorption { room.air.nepers_per_metre(BANDS[b]) } else { 0.0 }
    });
    let sphere = 4.0 / 3.0 * std::f64::consts::PI * rays.radius.powi(3);
    let each = 1.0 / rays.count as f64;
    let reach = c * rays.seconds;
    let empty = || vec![vec![[0.0f64; BANDS.len()]; bins]; receivers.len()];
    let histograms = (0..rays.count)
        .into_par_iter()
        .fold(empty, |mut acc, index| {
            stage.advance(1.0);
            let mut random = Random(rays.seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ index as u64);
            let mut origin = source;
            let mut direction = random.sphere();
            let mut energy = [each; BANDS.len()];
            let mut travelled = 0.0;
            let mut skip = usize::MAX;
            let mut order = 0;
            while travelled < reach && energy.iter().any(|e| *e > each * 1e-9) {
                let Some((length, hit)) = bvh.ray_hit(origin, direction, skip, 1e-9) else { break };
                if order > 0 {
                    for (r, centre) in receivers.iter().enumerate() {
                        let Some((enter, inside)) =
                            chord(origin, direction, length, *centre, rays.radius)
                        else {
                            continue;
                        };
                        let at = (travelled + enter + inside / 2.0) / c;
                        let bin = (at / rays.bin) as usize;
                        if bin >= bins {
                            continue;
                        }
                        for b in 0..BANDS.len() {
                            let fade = (-air[b] * (enter + inside / 2.0)).exp();
                            acc[r][bin][b] += energy[b] * fade * inside / (c * sphere);
                        }
                    }
                }
                for b in 0..BANDS.len() {
                    energy[b] *= (-air[b] * length).exp();
                }
                travelled += length;
                origin = add(origin, scale(direction, length));
                let surface = &room.surfaces[owner[hit]];
                energy.iter_mut().zip(surface.absorption).for_each(|(e, a)| *e *= 1.0 - a);
                let mut n = normals[hit];
                if dot(n, direction) > 0.0 {
                    n = scale(n, -1.0);
                }
                direction = if random.next() < surface.scattering {
                    random.lambert(n)
                } else {
                    sub(direction, scale(n, 2.0 * dot(direction, n)))
                };
                skip = hit;
                order += 1;
            }
            acc
        })
        .reduce(empty, |mut a, b| {
            for (x, y) in a.iter_mut().zip(b) {
                for (p, q) in x.iter_mut().zip(y) {
                    for k in 0..BANDS.len() {
                        p[k] += q[k];
                    }
                }
            }
            a
        });
    Ok(receivers
        .iter()
        .zip(histograms)
        .map(|(r, h)| {
            let d = dot(sub(*r, source), sub(*r, source)).sqrt();
            let direct_at = d / c;
            let mut energy: [Vec<f64>; BANDS.len()] =
                std::array::from_fn(|b| h.iter().map(|v| v[b]).collect());
            let bin = (direct_at / rays.bin) as usize;
            if bin < bins && d > 0.0 && visible(&bvh, source, *r) {
                for (b, e) in energy.iter_mut().enumerate() {
                    e[bin] += (-air[b] * d).exp() / (4.0 * std::f64::consts::PI * d * d * c);
                }
            }
            Response { bin: rays.bin, direct: direct_at, energy }
        })
        .collect())
}

fn visible(bvh: &Bvh, a: [f64; 3], b: [f64; 3]) -> bool {
    let d = sub(b, a);
    let length = dot(d, d).sqrt();
    bvh.ray(a, scale(d, 1.0 / length), usize::MAX, 1e-9) >= length
}

pub fn shoebox(size: [f64; 3], absorption: [f64; BANDS.len()], scattering: f64) -> Vec<Surface> {
    let [x, y, z] = size;
    let p = |i: usize| {
        [
            if i & 1 > 0 { x } else { 0.0 },
            if i & 2 > 0 { y } else { 0.0 },
            if i & 4 > 0 { z } else { 0.0 },
        ]
    };
    let quads =
        [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
    quads
        .iter()
        .map(|q| Surface {
            triangles: vec![[p(q[0]), p(q[1]), p(q[2])], [p(q[0]), p(q[2]), p(q[3])]],
            absorption,
            scattering,
        })
        .collect()
}
