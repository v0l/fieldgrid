pub type Point = [f64; 3];
pub type Triangle = [Point; 3];

struct Node {
    low: Point,
    high: Point,
    start: usize,
    count: usize,
    children: Option<(usize, usize)>,
}

pub struct Bvh {
    triangles: Vec<Triangle>,
    order: Vec<usize>,
    nodes: Vec<Node>,
}

fn sub(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Point, b: Point) -> Point {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn bounds_of(triangles: &[Triangle], items: &[usize]) -> (Point, Point) {
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for &i in items {
        for p in triangles[i] {
            for k in 0..3 {
                low[k] = low[k].min(p[k]);
                high[k] = high[k].max(p[k]);
            }
        }
    }
    (low, high)
}

impl Bvh {
    pub fn new(triangles: Vec<Triangle>) -> Bvh {
        let mut bvh = Bvh { order: (0..triangles.len()).collect(), triangles, nodes: Vec::new() };
        if !bvh.order.is_empty() {
            bvh.build(0, bvh.order.len());
        }
        bvh
    }

    pub fn indexed(points: &[Point], triangles: &[[usize; 3]]) -> Bvh {
        Bvh::new(triangles.iter().map(|t| t.map(|v| points[v])).collect())
    }

    pub fn triangles(&self) -> &[Triangle] {
        &self.triangles
    }

    pub fn bounds(&self) -> Option<(Point, Point)> {
        self.nodes.first().map(|n| (n.low, n.high))
    }

    fn build(&mut self, start: usize, end: usize) -> usize {
        let (low, high) = bounds_of(&self.triangles, &self.order[start..end]);
        let index = self.nodes.len();
        self.nodes.push(Node { low, high, start, count: end - start, children: None });
        if end - start <= 4 {
            return index;
        }
        let extent = sub(high, low);
        let axis = if extent[0] >= extent[1] && extent[0] >= extent[2] {
            0
        } else if extent[1] >= extent[2] {
            1
        } else {
            2
        };
        let triangles = &self.triangles;
        let centre = |i: usize| {
            let [a, b, c] = triangles[i];
            (a[axis] + b[axis] + c[axis]) / 3.0
        };
        let middle = (start + end) / 2;
        self.order[start..end]
            .select_nth_unstable_by(middle - start, |&a, &b| centre(a).total_cmp(&centre(b)));
        let left = self.build(start, middle);
        let right = self.build(middle, end);
        self.nodes[index].children = Some((left, right));
        index
    }

    fn slab(node: &Node, origin: Point, inverse: Point, limit: f64) -> bool {
        let mut near: f64 = 0.0;
        let mut far = limit;
        for k in 0..3 {
            let t1 = (node.low[k] - origin[k]) * inverse[k];
            let t2 = (node.high[k] - origin[k]) * inverse[k];
            let (a, b) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
            near = near.max(a);
            far = far.min(b);
            if near > far {
                return false;
            }
        }
        true
    }

    pub fn ray(&self, origin: Point, direction: Point, skip: usize, least: f64) -> f64 {
        let mut best = f64::INFINITY;
        if self.nodes.is_empty() {
            return best;
        }
        let inverse = direction.map(|d| 1.0 / d);
        let mut stack = vec![0];
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n];
            if !Self::slab(node, origin, inverse, best) {
                continue;
            }
            match node.children {
                Some((l, r)) => {
                    stack.push(l);
                    stack.push(r);
                }
                None => {
                    for &i in &self.order[node.start..node.start + node.count] {
                        if i == skip {
                            continue;
                        }
                        if let Some(t) = hit(self.triangles[i], origin, direction)
                            && t > least
                            && t < best
                        {
                            best = t;
                        }
                    }
                }
            }
        }
        best
    }
}

fn hit([a, b, c]: Triangle, origin: Point, direction: Point) -> Option<f64> {
    let (e1, e2) = (sub(b, a), sub(c, a));
    let p = cross(direction, e2);
    let det = dot(e1, p);
    if det.abs() < 1.0e-18 {
        return None;
    }
    let inv = 1.0 / det;
    let s = sub(origin, a);
    let u = dot(s, p) * inv;
    if !(-1.0e-9..=1.0 + 1.0e-9).contains(&u) {
        return None;
    }
    let q = cross(s, e1);
    let v = dot(direction, q) * inv;
    if v < -1.0e-9 || u + v > 1.0 + 1.0e-9 {
        return None;
    }
    Some(dot(e2, q) * inv)
}
