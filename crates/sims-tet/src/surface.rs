use crate::fixed::{HashMap, HashSet};
use crate::geom::*;

pub struct SurfaceMesh<S> {
    pub points: Vec<Point3>,
    pub triangles: Vec<[u32; 3]>,
    pub tags: Vec<usize>,
    pub uvs: Vec<Vec<(usize, Vector2)>>,
    pub faces: Vec<S>,
}

pub struct Welder {
    pub quantum: f64,
    pub index: HashMap<[i64; 3], u32>,
}

impl Welder {
    fn key(&self, p: Point3) -> [i64; 3] {
        [p.x, p.y, p.z].map(|v| (v / self.quantum).round() as i64)
    }
}

impl<S: Surface> SurfaceMesh<S> {
    pub fn welded(&mut self, welder: &mut Welder, p: Point3) -> u32 {
        let key = welder.key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let near = [key[0] + dx, key[1] + dy, key[2] + dz];
                    if let Some(&id) = welder.index.get(&near)
                        && (self.points[id as usize] - p).magnitude() <= welder.quantum
                    {
                        return id;
                    }
                }
            }
        }
        let id = self.fresh(p);
        welder.index.insert(key, id);
        id
    }

    pub fn fresh(&mut self, p: Point3) -> u32 {
        self.points.push(p);
        self.uvs.push(Vec::new());
        (self.points.len() - 1) as u32
    }

    pub fn note_uv(&mut self, id: u32, face: usize, uv: Vector2) {
        let list = &mut self.uvs[id as usize];
        if !list.iter().any(|(f, _)| *f == face) {
            list.push((face, uv));
        }
    }

    pub fn uv_on(&self, id: u32, face: usize) -> Option<Vector2> {
        self.uvs[id as usize].iter().find(|(f, _)| *f == face).map(|(_, uv)| *uv)
    }

    pub fn edges(&self) -> HashMap<[u32; 2], Vec<usize>> {
        let mut edges: HashMap<[u32; 2], Vec<usize>> = HashMap::default();
        for (t, tri) in self.triangles.iter().enumerate() {
            for k in 0..3 {
                let mut e = [tri[k], tri[(k + 1) % 3]];
                e.sort_unstable();
                edges.entry(e).or_default().push(t);
            }
        }
        edges
    }

    pub fn open_edges(&self) -> usize {
        self.edges().values().filter(|t| t.len() != 2).count()
    }

    pub fn delaunay_flips(&mut self) -> usize {
        let mut edges = self.edges();
        let mut queue: Vec<[u32; 2]> = edges.keys().copied().collect();
        queue.sort_unstable();
        let mut flips = 0;
        let limit = self.triangles.len() * 20;
        while let Some(edge) = queue.pop() {
            if flips > limit {
                break;
            }
            let Some(owners) = edges.get(&edge) else {
                continue;
            };
            let [t1, t2] = owners[..] else { continue };
            if self.tags[t1] != self.tags[t2] {
                continue;
            }
            let [a, b] = edge;
            let apex = |t: usize| {
                *self.triangles[t]
                    .iter()
                    .find(|&&v| v != a && v != b)
                    .expect("a triangle has a third vertex")
            };
            let (c, d) = (apex(t1), apex(t2));
            if c == d {
                continue;
            }
            let mut diagonal = [c, d];
            diagonal.sort_unstable();
            if edges.contains_key(&diagonal) {
                continue;
            }
            let p = |v: u32| self.points[v as usize];
            let angle = |at: u32| {
                let (u, w) = (p(a) - p(at), p(b) - p(at));
                (u.dot(w) / (u.magnitude() * w.magnitude()).max(1.0e-300)).clamp(-1.0, 1.0).acos()
            };
            if angle(c) + angle(d) <= std::f64::consts::PI + 1.0e-9 {
                continue;
            }
            let tri = self.triangles[t1];
            let k = (0..3).find(|&k| tri[k] == c).expect("apex");
            let (first, second) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
            let new1 = [c, first, d];
            let new2 = [d, second, c];
            let normal = |t: [u32; 3]| (p(t[1]) - p(t[0])).cross(p(t[2]) - p(t[0]));
            let old = normal(tri);
            if normal(new1).dot(old) <= 0.0 || normal(new2).dot(old) <= 0.0 {
                continue;
            }
            for (t, tri) in [(t1, self.triangles[t1]), (t2, self.triangles[t2])] {
                for k in 0..3 {
                    let mut e = [tri[k], tri[(k + 1) % 3]];
                    e.sort_unstable();
                    if let Some(list) = edges.get_mut(&e) {
                        list.retain(|&x| x != t);
                        if list.is_empty() {
                            edges.remove(&e);
                        }
                    }
                }
            }
            self.triangles[t1] = new1;
            self.triangles[t2] = new2;
            for (t, tri) in [(t1, new1), (t2, new2)] {
                for k in 0..3 {
                    let mut e = [tri[k], tri[(k + 1) % 3]];
                    e.sort_unstable();
                    edges.entry(e).or_default().push(t);
                    if e != diagonal {
                        queue.push(e);
                    }
                }
            }
            flips += 1;
        }
        flips
    }

    pub fn split_edge(&mut self, edge: [u32; 2], owners: &[usize]) -> Point3 {
        let [a, b] = edge;
        let same = owners.iter().map(|&t| self.tags[t]).collect::<HashSet<_>>();
        let (pa, pb) = (self.points[a as usize], self.points[b as usize]);
        let mut at = pa.midpoint(pb);
        let mut uv_note = None;
        if same.len() == 1 {
            let face = self.tags[owners[0]];
            if let (Some(ua), Some(ub)) = (self.uv_on(a, face), self.uv_on(b, face)) {
                let uv = (ua + ub) / 2.0;
                if let Some(on) = self.faces[face].at(uv)
                    && (on - at).magnitude() < (pb - pa).magnitude() * 0.25
                {
                    at = on;
                    uv_note = Some((face, uv));
                }
            }
        }
        let m = self.fresh(at);
        if let Some((face, uv)) = uv_note {
            self.note_uv(m, face, uv);
        }
        for &t in owners {
            let tri = self.triangles[t];
            let k = (0..3)
                .find(|&k| {
                    let e = [tri[k], tri[(k + 1) % 3]];
                    e == [a, b] || e == [b, a]
                })
                .expect("the owner holds the edge");
            let (p, q, r) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
            self.triangles[t] = [p, m, r];
            self.triangles.push([m, q, r]);
            self.tags.push(self.tags[t]);
        }
        at
    }
}
