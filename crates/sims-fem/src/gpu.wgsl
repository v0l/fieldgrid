struct Range {
    first: u32,
    count: u32,
    nodes: u32,
    pad: u32,
}

@group(0) @binding(0) var<uniform> range: Range;
@group(0) @binding(1) var<storage, read> index_a: array<u32>;
@group(0) @binding(2) var<storage, read> index_b: array<u32>;
@group(0) @binding(3) var<storage, read> matrix: array<f64>;
@group(0) @binding(4) var<storage, read> source: array<f64>;
@group(0) @binding(5) var<storage, read_write> output: array<f64>;
@group(0) @binding(6) var<storage, read> extra: array<f64>;
@group(0) @binding(7) var<storage, read> scalars: array<f64>;

fn flat(g: vec3<u32>) -> u32 {
    return g.x + g.y * 65535u * 256u;
}

fn block_times(base: u32, v: vec3<f64>) -> vec3<f64> {
    return vec3<f64>(
        matrix[base] * v.x + matrix[base + 1u] * v.y + matrix[base + 2u] * v.z,
        matrix[base + 3u] * v.x + matrix[base + 4u] * v.y + matrix[base + 5u] * v.z,
        matrix[base + 6u] * v.x + matrix[base + 7u] * v.y + matrix[base + 8u] * v.z,
    );
}

fn source_at(node: u32) -> vec3<f64> {
    return vec3<f64>(source[3u * node], source[3u * node + 1u], source[3u * node + 2u]);
}

@compute @workgroup_size(256, 1, 1)
fn multiply(@builtin(global_invocation_id) g: vec3<u32>) {
    let local = flat(g);
    if local >= range.count { return; }
    var sum = vec3<f64>(0.0, 0.0, 0.0);
    let slice = local / 32u;
    let lane = local % 32u;
    let start = index_a[2u * slice];
    let width = index_a[2u * slice + 1u];
    for (var k = 0u; k < width; k++) {
        let e = start + 32u * k;
        let v = source_at(index_b[e + lane]);
        let m = 9u * e + lane;
        sum += vec3<f64>(
            matrix[m] * v.x + matrix[m + 32u] * v.y + matrix[m + 64u] * v.z,
            matrix[m + 96u] * v.x + matrix[m + 128u] * v.y + matrix[m + 160u] * v.z,
            matrix[m + 192u] * v.x + matrix[m + 224u] * v.y + matrix[m + 256u] * v.z,
        );
    }
    let row = range.first + local;
    output[3u * row] = sum.x;
    output[3u * row + 1u] = sum.y;
    output[3u * row + 2u] = sum.z;
}

@compute @workgroup_size(256, 1, 1)
fn project(@builtin(global_invocation_id) g: vec3<u32>) {
    let node = flat(g);
    if node >= range.nodes { return; }
    let v = vec3<f64>(output[3u * node], output[3u * node + 1u], output[3u * node + 2u]);
    let p = block_times(9u * node, v);
    output[3u * node] = p.x;
    output[3u * node + 1u] = p.y;
    output[3u * node + 2u] = p.z;
}

@compute @workgroup_size(256, 1, 1)
fn sweep(@builtin(global_invocation_id) g: vec3<u32>) {
    let node = flat(g);
    if node >= range.nodes { return; }
    let p = block_times(9u * node, source_at(node));
    output[3u * node] = extra[3u * node] - p.x;
    output[3u * node + 1u] = extra[3u * node + 1u] - p.y;
    output[3u * node + 2u] = extra[3u * node + 2u] - p.z;
}

@compute @workgroup_size(256, 1, 1)
fn jacobi(@builtin(global_invocation_id) g: vec3<u32>) {
    let node = flat(g);
    if node >= range.nodes { return; }
    let w = scalars[range.first];
    let v = block_times(9u * node, source_at(node));
    output[3u * node] += w * v.x;
    output[3u * node + 1u] += w * v.y;
    output[3u * node + 2u] += w * v.z;
}

@compute @workgroup_size(256, 1, 1)
fn clusters(@builtin(global_invocation_id) g: vec3<u32>) {
    let row = flat(g);
    if row >= range.count { return; }
    let first_node = index_a[4u * row];
    let first_entry = index_a[4u * row + 1u];
    let size = index_a[4u * row + 2u];
    let own = index_a[4u * row + 3u];
    var sum: f64 = 0.0;
    for (var c = 0u; c < size; c++) {
        let node = index_b[first_node + c / 3u];
        sum += matrix[first_entry + c] * source[3u * node + c % 3u];
    }
    let node = index_b[first_node + own / 3u];
    output[3u * node + own % 3u] += scalars[range.first] * sum;
}

@compute @workgroup_size(256, 1, 1)
fn gather(@builtin(global_invocation_id) g: vec3<u32>) {
    let coarse = flat(g);
    if coarse >= range.count { return; }
    var sum = vec3<f64>(0.0, 0.0, 0.0);
    for (var k = index_a[coarse]; k < index_a[coarse + 1u]; k++) {
        sum += matrix[k] * source_at(index_b[k]);
    }
    output[3u * coarse] = sum.x;
    output[3u * coarse + 1u] = sum.y;
    output[3u * coarse + 2u] = sum.z;
}

@compute @workgroup_size(256, 1, 1)
fn prolong(@builtin(global_invocation_id) g: vec3<u32>) {
    let node = flat(g);
    if node >= range.nodes { return; }
    for (var k = 0u; k < 2u; k++) {
        let w = matrix[2u * node + k];
        if w != 0.0 {
            let v = source_at(index_a[2u * node + k]);
            output[3u * node] += w * v.x;
            output[3u * node + 1u] += w * v.y;
            output[3u * node + 2u] += w * v.z;
        }
    }
}

@compute @workgroup_size(256, 1, 1)
fn add_scaled(@builtin(global_invocation_id) g: vec3<u32>) {
    let i = flat(g);
    if i >= 3u * range.nodes { return; }
    output[i] += scalars[range.first] * source[i];
}

@compute @workgroup_size(256, 1, 1)
fn scale_add(@builtin(global_invocation_id) g: vec3<u32>) {
    let i = flat(g);
    if i >= 3u * range.nodes { return; }
    output[i] = source[i] + scalars[range.first] * output[i];
}

@compute @workgroup_size(256, 1, 1)
fn aggregate(@builtin(global_invocation_id) g: vec3<u32>) {
    let t = flat(g);
    if t >= 6u * range.count { return; }
    let group = t / 6u;
    let j = t % 6u;
    var sum: f64 = 0.0;
    for (var k = index_a[group]; k < index_a[group + 1u]; k++) {
        let v = index_b[k];
        for (var i = 0u; i < 3u; i++) {
            sum += matrix[18u * v + 6u * i + j] * source[3u * v + i];
        }
    }
    output[t] = sum;
}

@compute @workgroup_size(256, 1, 1)
fn spread(@builtin(global_invocation_id) g: vec3<u32>) {
    let v = flat(g);
    if v >= range.nodes { return; }
    let group = index_a[v];
    for (var i = 0u; i < 3u; i++) {
        var sum: f64 = 0.0;
        for (var j = 0u; j < 6u; j++) {
            sum += matrix[18u * v + 6u * i + j] * source[6u * group + j];
        }
        output[3u * v + i] += sum;
    }
}

var<workgroup> partial: array<f64, 256>;

@compute @workgroup_size(256, 1, 1)
fn dot(@builtin(local_invocation_id) l: vec3<u32>, @builtin(workgroup_id) w: vec3<u32>) {
    let n = 3u * range.nodes;
    var sum: f64 = 0.0;
    var i = w.x * 256u + l.x;
    while i < n {
        sum += source[i] * extra[i];
        i += range.count * 256u;
    }
    partial[l.x] = sum;
    workgroupBarrier();
    for (var k = 128u; k > 0u; k /= 2u) {
        if l.x < k { partial[l.x] += partial[l.x + k]; }
        workgroupBarrier();
    }
    if l.x == 0u { output[w.x] = partial[0]; }
}
