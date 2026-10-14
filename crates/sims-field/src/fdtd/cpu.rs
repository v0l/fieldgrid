use super::run::{Packed, pack};
use super::{Extras, MU0, Pulse, Record, Sim, psi_idx};
use rayon::prelude::*;

#[derive(Clone, Copy)]
struct Shared(*mut f32);

unsafe impl Send for Shared {}
unsafe impl Sync for Shared {}

impl Shared {
    fn get(self) -> *mut f32 {
        self.0
    }
}

struct State<'a> {
    sim: &'a Sim,
    p: Packed,
    e: Vec<f32>,
    h: Vec<f32>,
    psi: Vec<f32>,
    inductor: Vec<f32>,
    sheets: Vec<f32>,
    current: Vec<f32>,
    debye_edges: Vec<f32>,
    debye_state: Vec<f32>,
    series: Vec<f32>,
    plane: Vec<f32>,
    ntff: Vec<f32>,
    patch_idx: Vec<u32>,
    omega: Vec<f32>,
    plane_k: usize,
    w0: f32,
}

fn shift(p: [usize; 3], a: usize, d: isize) -> [usize; 3] {
    let mut q = p;
    q[a] = (q[a] as isize + d) as usize;
    q
}

impl State<'_> {
    fn psi_base(&self, e: usize, c: usize, a: usize) -> usize {
        self.p.offsets[e * 9 + c * 3 + a] as usize
    }

    fn update_h(&mut self) {
        let sim = self.sim;
        let n = self.p.n;
        let nn = self.p.nn;
        let pml = sim.pml;
        let k_mu = (sim.dt / MU0) as f32;
        let e = &self.e;
        let h = Shared(self.h.as_mut_ptr());
        let psi = Shared(self.psi.as_mut_ptr());
        let bases: [[usize; 3]; 3] =
            std::array::from_fn(|c| std::array::from_fn(|a| self.psi_base(1, c, a)));
        let ax = &sim.ax;
        (0..n[0]).into_par_iter().for_each(|i| {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let p = [i, j, k];
                    let id = (i * n[1] + j) * n[2] + k;
                    for (c, base) in bases.iter().enumerate() {
                        let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                        if p[u] + 1 >= n[u] || p[v] + 1 >= n[v] {
                            continue;
                        }
                        let up = shift(p, u, 1);
                        let vp = shift(p, v, 1);
                        let at = |q: [usize; 3]| (q[0] * n[1] + q[1]) * n[2] + q[2];
                        let ev = e[v * nn + at(up)] - e[v * nn + id];
                        let eu = e[u * nn + at(vp)] - e[u * nn + id];
                        let (pu, pv) = (p[u], p[v]);
                        let mut t1 = ev * ax[u].inv_d[pu];
                        let mut t2 = eu * ax[v].inv_d[pv];
                        if pml > 0 {
                            let su = ax[u].slot_h[pu];
                            if su >= 0 {
                                let q = base[u] + psi_idx(n, u, pml, p, su as usize);
                                unsafe {
                                    let r = psi.get().add(q);
                                    *r =
                                        ax[u].bh[pu] * *r + ax[u].ch[pu] * ev * ax[u].raw_inv_d[pu];
                                    t1 += *r;
                                }
                            }
                            let sv = ax[v].slot_h[pv];
                            if sv >= 0 {
                                let q = base[v] + psi_idx(n, v, pml, p, sv as usize);
                                unsafe {
                                    let r = psi.get().add(q);
                                    *r =
                                        ax[v].bh[pv] * *r + ax[v].ch[pv] * eu * ax[v].raw_inv_d[pv];
                                    t2 += *r;
                                }
                            }
                        }
                        unsafe {
                            *h.get().add(c * nn + id) -= k_mu * (t1 - t2);
                        }
                    }
                }
            }
        });
    }

    fn update_e(&mut self) {
        let sim = self.sim;
        let n = self.p.n;
        let nn = self.p.nn;
        let pml = sim.pml;
        let h = &self.h;
        let e = Shared(self.e.as_mut_ptr());
        let psi = Shared(self.psi.as_mut_ptr());
        let bases: [[usize; 3]; 3] =
            std::array::from_fn(|c| std::array::from_fn(|a| self.psi_base(0, c, a)));
        let ax = &sim.ax;
        (0..n[0]).into_par_iter().for_each(|i| {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let p = [i, j, k];
                    let id = (i * n[1] + j) * n[2] + k;
                    for (c, base) in bases.iter().enumerate() {
                        let cb = sim.cb[c][id];
                        if cb == 0.0 {
                            continue;
                        }
                        let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                        let at = |q: [usize; 3]| (q[0] * n[1] + q[1]) * n[2] + q[2];
                        let hv = h[v * nn + id] - h[v * nn + at(shift(p, u, -1))];
                        let hu = h[u * nn + id] - h[u * nn + at(shift(p, v, -1))];
                        let (pu, pv) = (p[u], p[v]);
                        let mut t1 = hv * ax[u].inv_dd[pu];
                        let mut t2 = hu * ax[v].inv_dd[pv];
                        if pml > 0 {
                            let su = ax[u].slot_e[pu];
                            if su >= 0 {
                                let q = base[u] + psi_idx(n, u, pml, p, su as usize);
                                unsafe {
                                    let r = psi.get().add(q);
                                    *r = ax[u].be[pu] * *r
                                        + ax[u].ce[pu] * hv * ax[u].raw_inv_dd[pu];
                                    t1 += *r;
                                }
                            }
                            let sv = ax[v].slot_e[pv];
                            if sv >= 0 {
                                let q = base[v] + psi_idx(n, v, pml, p, sv as usize);
                                unsafe {
                                    let r = psi.get().add(q);
                                    *r = ax[v].be[pv] * *r
                                        + ax[v].ce[pv] * hu * ax[v].raw_inv_dd[pv];
                                    t2 += *r;
                                }
                            }
                        }
                        unsafe {
                            let f = e.get().add(c * nn + id);
                            *f = sim.ca[c][id] * *f + cb * (t1 - t2);
                        }
                    }
                }
            }
        });
    }

    fn sheet_before(&mut self, s: usize, old: f32, ha: f32, hb: f32) {
        let b = s * 16;
        let sh = &mut self.sheets;
        sh[b + 12] = old;
        let aa = sh[b + 7] + ha * ha;
        let bb = sh[b + 8] + hb * hb;
        let ab = sh[b + 9] + ha * hb;
        sh[b + 7] = aa;
        sh[b + 8] = bb;
        sh[b + 9] = ab;
        let sum = aa + bb;
        if sum <= 1e-30 {
            return;
        }
        let mut scale = 1.0;
        if sh[b + 6] != 0.0 {
            scale = sh[b + 10] + sh[b + 11];
            let total = aa + bb - 2.0 * ab;
            if total > 1e-30 {
                scale = (sum / total).clamp(0.5, 4.0);
            }
        }
        sh[b + 10] = scale * aa / sum;
        sh[b + 11] = scale * bb / sum;
    }

    fn sheet_after(&mut self, s: usize) {
        let b = s * 16;
        let branch = &self.sim.sheet_table;
        let sh = &mut self.sheets;
        let f = sh[b].to_bits() as usize;
        let (len, c, r) = (sh[b + 3], sh[b + 4], sh[b + 5]);
        let (ma, mb, old, i_old) = (sh[b + 10], sh[b + 11], sh[b + 12], sh[b + 13]);
        let nr = branch[0] as usize;
        let nc = branch[1] as usize;
        let stride = 2 + nr + 2 * nc;
        let faces = 4 + 2 * nr + 4 * nc;
        let fa = faces + sh[b + 14].to_bits() as usize * stride;
        let fb = faces + sh[b + 15].to_bits() as usize * stride;
        let base = s * (nr + 2 * nc);
        let cur = &mut self.current;
        let mut hist = (ma * branch[fa + 1] + mb * branch[fb + 1]) * i_old;
        for k in 0..nr {
            hist -= (ma * branch[fa + 2 + k] + mb * branch[fb + 2 + k]) * cur[base + k];
        }
        for k in 0..nc {
            let o = 2 + nr + 2 * k;
            let cr = ma * branch[fa + o] + mb * branch[fb + o];
            let ci = ma * branch[fa + o + 1] + mb * branch[fb + o + 1];
            let u = base + nr + 2 * k;
            hist -= cr * cur[u] - ci * cur[u + 1];
        }
        let a = r * (ma * branch[fa] + mb * branch[fb]);
        let e = &mut self.e;
        let ibar = (0.5 * len * (e[f] + old) + r * hist) / (a + 0.5 * len * c);
        e[f] -= c * ibar;
        sh[b + 13] = 2.0 * ibar - i_old;
        for k in 0..nr {
            let p = 4 + 2 * k;
            cur[base + k] = branch[p] * cur[base + k] + branch[p + 1] * ibar;
        }
        for k in 0..nc {
            let p = 4 + 2 * nr + 4 * k;
            let u = base + nr + 2 * k;
            let (ur, ui) = (cur[u], cur[u + 1]);
            cur[u] = branch[p] * ur - branch[p + 1] * ui + branch[p + 2] * ibar;
            cur[u + 1] = branch[p] * ui + branch[p + 1] * ur + branch[p + 3] * ibar;
        }
    }

    fn debye_after(&mut self, s: usize) {
        let f = self.debye_edges[s * 3].to_bits() as usize;
        let dd = self.debye_edges[s * 3 + 1];
        let old = self.debye_edges[s * 3 + 2];
        let nn = self.p.nn;
        let cb = self.sim.cb[f / nn][f % nn];
        let table = &self.p.debye_table;
        let np = table[0] as usize;
        let count = self.p.debye_count;
        let mut s0 = 0.0;
        for k in 0..np {
            s0 += 0.5 * (1.0 + table[1 + 2 * k]) * self.debye_state[k * count + s];
        }
        let next = self.e[f] - cb * s0;
        self.e[f] = next;
        let de = next - old;
        for k in 0..np {
            let q = k * count + s;
            self.debye_state[q] =
                table[1 + 2 * k] * self.debye_state[q] + dd * table[2 + 2 * k] * de;
        }
    }

    fn pulse(&self, t: f32) -> f32 {
        let x = (t - self.p.t0 as f32) / self.p.tau as f32;
        (-x * x).exp() * (self.w0 * (t - self.p.t0 as f32)).cos()
    }

    fn step(&mut self, n: usize) {
        let nn = self.p.nn;
        self.update_h();
        for w in self.p.wire.as_chunks::<9>().0 {
            let e = |k: usize| f32::from_bits(w[k + 1]) * self.e[w[k] as usize];
            self.h[w[0] as usize] += e(1) + e(3) + e(5) + e(7);
        }
        for s in 0..self.p.sheet_count {
            let b = s * 16;
            let f = self.sheets[b].to_bits() as usize;
            let ha = self.h[self.sheets[b + 1].to_bits() as usize];
            let hb = self.h[self.sheets[b + 2].to_bits() as usize];
            self.sheet_before(s, self.e[f], ha, hb);
        }
        for s in 0..self.p.debye_count {
            self.debye_edges[s * 3 + 2] = self.e[self.debye_edges[s * 3].to_bits() as usize];
        }
        self.update_e();
        for s in 0..self.p.sheet_count {
            self.sheet_after(s);
        }
        for s in 0..self.p.debye_count {
            self.debye_after(s);
        }
        let dt = self.sim.dt as f32;
        let t = (n as f32 + 0.5) * dt;
        let pulse = self.pulse(t);
        for q in self.p.source.as_chunks::<3>().0 {
            let (id, src, comp) = (q[0] as usize, q[1], q[2] as usize);
            self.e[comp * nn + id] -= src * pulse;
        }
        for q in self.inductor.as_chunks_mut::<5>().0 {
            let f = q[0] as usize * nn + q[1] as usize;
            let x = self.e[f] - q[2] * q[4];
            self.e[f] = x;
            q[4] += q[3] * x;
        }
        self.probe(n);
        self.field_dft(n);
        self.far(n);
    }

    fn probe(&mut self, n: usize) {
        let nn = self.p.nn;
        let ports = self.p.ports;
        for k in 0..ports {
            let d = &self.p.probe_def[k * 8..k * 8 + 8];
            let comp = d[0] as usize;
            let (first, count) = (d[1] as usize, d[2] as usize);
            let (lf, lc) = (d[3] as usize, d[4] as usize);
            let mut v = 0.0f64;
            for e in first..first + count {
                let q = e * 2;
                v -= (self.e[comp * nn + self.p.port_edges[q] as usize] * self.p.port_edges[q + 1])
                    as f64;
            }
            let mut cur = 0.0f64;
            for c in lf..lf + lc {
                let q = c * 3;
                let hc = self.p.loops[q + 1] as usize - 3;
                cur += (self.h[hc * nn + self.p.loops[q] as usize] * self.p.loops[q + 2]) as f64;
            }
            if n < self.p.cap {
                self.series[(n * ports + k) * 2] = v as f32;
                self.series[(n * ports + k) * 2 + 1] = cur as f32;
            }
        }
    }

    fn field_dft(&mut self, n: usize) {
        if self.p.plane_len == 0 {
            return;
        }
        let [n0, n1, n2] = self.p.n;
        let nn = self.p.nn;
        let nf = self.p.nf;
        let dt = self.sim.dt as f32;
        let te = (n as f32 + 1.0) * dt;
        let th = (n as f32 + 0.5) * dt;
        let (e, h, omega, plane_k) = (&self.e, &self.h, &self.omega, self.plane_k);
        self.plane.par_chunks_mut(nf * 10).enumerate().take(n0 * n1).for_each(|(q, out)| {
            let (i, j) = (q / n1, q % n1);
            let id = (i * n1 + j) * n2 + plane_k;
            let f = [e[id], e[nn + id], e[2 * nn + id], h[id], h[nn + id]];
            for (k, w) in omega.iter().enumerate().take(nf) {
                let ce = ((w * te).cos() * dt, -(w * te).sin() * dt);
                let ch = ((w * th).cos() * dt, -(w * th).sin() * dt);
                for (c, v) in f.iter().enumerate() {
                    let ph = if c >= 3 { ch } else { ce };
                    out[k * 10 + 2 * c] += v * ph.0;
                    out[k * 10 + 2 * c + 1] += v * ph.1;
                }
            }
        });
    }

    fn far(&mut self, n: usize) {
        if self.p.patches == 0 || self.p.nf == 0 {
            return;
        }
        let nn = self.p.nn;
        let nf = self.p.nf;
        let dt = self.sim.dt as f32;
        let te = (n as f32 + 1.0) * dt;
        let th = (n as f32 + 0.5) * dt;
        let (e, h, idx, omega) = (&self.e, &self.h, &self.patch_idx, &self.omega);
        self.ntff.par_chunks_mut(nf * 12).enumerate().for_each(|(g, out)| {
            let b = g * 13;
            let meta = idx[b];
            let axis = (meta & 3) as usize;
            let side = if meta & 4 != 0 { 1.0 } else { -1.0 };
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            let at = |k: usize| idx[b + k] as usize;
            let mut ef = [0f32; 3];
            let mut hf = [0f32; 3];
            ef[u] = 0.5 * (e[u * nn + at(1)] + e[u * nn + at(2)]);
            ef[v] = 0.5 * (e[v * nn + at(3)] + e[v * nn + at(4)]);
            hf[u] = 0.25
                * (h[u * nn + at(5)] + h[u * nn + at(6)] + h[u * nn + at(7)] + h[u * nn + at(8)]);
            hf[v] = 0.25
                * (h[v * nn + at(9)]
                    + h[v * nn + at(10)]
                    + h[v * nn + at(11)]
                    + h[v * nn + at(12)]);
            let mut nrm = [0f32; 3];
            nrm[axis] = side;
            let jv = [
                nrm[1] * hf[2] - nrm[2] * hf[1],
                nrm[2] * hf[0] - nrm[0] * hf[2],
                nrm[0] * hf[1] - nrm[1] * hf[0],
            ];
            let mv = [
                ef[1] * nrm[2] - ef[2] * nrm[1],
                ef[2] * nrm[0] - ef[0] * nrm[2],
                ef[0] * nrm[1] - ef[1] * nrm[0],
            ];
            for (k, wf) in omega.iter().enumerate().take(nf) {
                let pe = ((-wf * te).cos() * dt, (-wf * te).sin() * dt);
                let ph = ((-wf * th).cos() * dt, (-wf * th).sin() * dt);
                let o = k * 12;
                for c in 0..3 {
                    out[o + 2 * c] += ph.0 * jv[c];
                    out[o + 2 * c + 1] += ph.1 * jv[c];
                    out[o + 6 + 2 * c] += pe.0 * mv[c];
                    out[o + 6 + 2 * c + 1] += pe.1 * mv[c];
                }
            }
        });
    }

    fn energy(&self) -> f64 {
        let e: f64 = self.e.par_iter().map(|x| (*x as f64).powi(2)).sum();
        let h: f64 = self.h.par_iter().map(|x| (*x as f64 * 376.73).powi(2)).sum();
        e + h
    }
}

pub fn run_cpu(
    sim: &Sim,
    driven: usize,
    pulse: &Pulse,
    extras: &Extras,
    progress: &mut dyn FnMut(usize, f64) -> bool,
) -> Result<Record, String> {
    let p = pack(sim, driven, pulse, extras)?;
    let nn = p.nn;
    let omega: Vec<f32> =
        extras.freqs.iter().take(4).map(|f| (2.0 * std::f64::consts::PI * f) as f32).collect();
    let mut st = State {
        sim,
        e: vec![0.0; 3 * nn],
        h: vec![0.0; 3 * nn],
        psi: vec![0.0; p.psi_total],
        inductor: p.inductor.clone(),
        sheets: p.sheets.clone(),
        current: vec![0.0; p.sheet_count * p.states],
        debye_edges: p.debye_edges.clone(),
        debye_state: vec![0.0; p.debye_count * sim.debye.poles.len()],
        series: vec![0.0; p.cap * p.ports * 2],
        plane: vec![0.0; p.plane_len],
        ntff: vec![0.0; p.patches * p.nf * 12],
        patch_idx: extras.ntff.clone().unwrap_or_default(),
        omega,
        plane_k: extras.plane_k.unwrap_or(0),
        w0: (2.0 * std::f64::consts::PI * pulse.f0) as f32,
        p,
    };
    let chunk = if extras.check_every == 0 { 1000 } else { extras.check_every };
    let pulse_end = (2.0 * st.p.t0 / sim.dt) as usize;
    let (mut steps, mut peak, mut decay) = (0usize, 0f64, 0f64);
    let stage = sims_progress::stage("FDTD time steps", Some(extras.max_steps), "steps");
    while steps < extras.max_steps {
        let k = chunk.min(extras.max_steps - steps);
        for s in 0..k {
            st.step(steps + s);
        }
        steps += k;
        let e = st.energy();
        if !e.is_finite() {
            return Err("the FDTD run went unstable".into());
        }
        peak = peak.max(e);
        decay = if e > 0.0 && peak > 0.0 { 10.0 * (peak / e).log10() } else { 0.0 };
        let goal = extras.decay_db;
        stage.progress(steps as f64, || format!("{decay:.1} of {goal} dB down"));
        if !progress(steps, decay) {
            return Err("cancelled".into());
        }
        if steps > pulse_end.max(extras.min_steps) && decay >= extras.decay_db {
            break;
        }
    }
    let ports = st.p.ports;
    st.series.truncate(steps * ports * 2);
    Ok(Record { steps, ports, series: st.series, decay_db: decay, plane: st.plane, ntff: st.ntff })
}
