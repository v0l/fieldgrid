use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type Node = Option<usize>;

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ElemKind {
    Resistor {
        r: f64,
    },
    Capacitor {
        c: f64,
    },
    Inductor {
        l: f64,
    },
    Diode {
        is: f64,
        n: f64,
        bv: Option<f64>,
        cj: f64,
    },
    Bjt {
        pnp: bool,
        is: f64,
        bf: f64,
        br: f64,
        va: Option<f64>,
        cje: f64,
        cjc: f64,
    },
    Mosfet {
        p: bool,
        vth: f64,
        k: f64,
        lambda: f64,
        cgs: f64,
        cgd: f64,
    },
    Opamp {
        gain: f64,
        gbw: f64,
        rout: f64,
        vos: f64,
        en: f64,
        fc: f64,
        headroom_high: f64,
        headroom_low: f64,
        slew: f64,
        ib: f64,
        psrr: f64,
    },
    Ldo {
        vout: f64,
        dropout: f64,
        rout: f64,
        iq: f64,
    },
    Vcvs {
        gain: f64,
    },
    Vccs {
        gm: f64,
    },
    Varistor {
        v1ma: f64,
        alpha: f64,
        c: f64,
    },
}

#[derive(Clone, Debug, Serialize)]
pub struct Elem {
    pub name: String,
    pub kind: ElemKind,
    pub nodes: Vec<Node>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Wave {
    Dc,
    Pulse { low: f64, high: f64, delay: f64, rise: f64, fall: f64, width: f64, period: f64 },
    Sine { offset: f64, amplitude: f64, frequency: f64, delay: f64 },
    Steps { points: Vec<(f64, f64)> },
    Surge { offset: f64, peak: f64, delay: f64, shape: Surge },
    Burst { offset: f64, peak: f64, delay: f64, rate: f64, length: f64, period: f64 },
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Surge {
    Open,
    Short,
}

struct Heidler {
    rise: f64,
    fall: f64,
    order: f64,
}

const OPEN: Heidler = Heidler { rise: 0.356e-6, fall: 65.845e-6, order: 1.852 };
const SHORT: Heidler = Heidler { rise: 47.52e-6, fall: 4.296e-6, order: 2.741 };
const FAST: Heidler = Heidler { rise: 3.5e-9, fall: 50e-9, order: 1.8 };

impl Heidler {
    fn raw(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let x = (t / self.rise).powf(self.order);
        x / (1.0 + x) * (-t / self.fall).exp()
    }

    fn peak(&self) -> f64 {
        let (mut lo, mut hi) = (0.0, 20.0 * self.rise.max(self.fall));
        for _ in 0..200 {
            let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
            if self.raw(a) < self.raw(b) { lo = a } else { hi = b }
        }
        self.raw(0.5 * (lo + hi))
    }

    fn at(&self, t: f64) -> f64 {
        self.raw(t) / self.peak()
    }

    fn span(&self) -> f64 {
        self.rise.min(self.fall)
    }
}

fn burst_pulse(u: f64, rate: f64, length: f64) -> Option<f64> {
    if u < 0.0 || u >= length {
        return None;
    }
    Some(u - (u * rate).floor() / rate)
}

impl Surge {
    fn heidler(self) -> &'static Heidler {
        match self {
            Surge::Open => &OPEN,
            Surge::Short => &SHORT,
        }
    }
}

fn marks(start: f64, first: f64, last: f64, stop: f64, out: &mut Vec<f64>) {
    out.push(start);
    let mut e = first;
    while e <= last && start + e <= stop {
        out.push(start + e);
        e *= std::f64::consts::SQRT_2;
    }
}

impl Wave {
    pub fn at(&self, dc: f64, t: f64) -> f64 {
        match self {
            Wave::Dc => dc,
            Wave::Pulse { low, high, delay, rise, fall, width, period } => {
                if t < *delay {
                    return *low;
                }
                let mut u = t - delay;
                if *period > 0.0 {
                    u %= period;
                }
                if u < *rise {
                    low + (high - low) * u / rise.max(1e-30)
                } else if u < rise + width {
                    *high
                } else if u < rise + width + fall {
                    high - (high - low) * (u - rise - width) / fall.max(1e-30)
                } else {
                    *low
                }
            }
            Wave::Sine { offset, amplitude, frequency, delay } => {
                if t < *delay {
                    *offset
                } else {
                    offset + amplitude * (std::f64::consts::TAU * frequency * (t - delay)).sin()
                }
            }
            Wave::Steps { points } => {
                let Some(first) = points.first() else { return dc };
                if t <= first.0 {
                    return first.1;
                }
                for w in points.windows(2) {
                    if t <= w[1].0 {
                        let span = (w[1].0 - w[0].0).max(1e-30);
                        return w[0].1 + (w[1].1 - w[0].1) * (t - w[0].0) / span;
                    }
                }
                points.last().map_or(dc, |p| p.1)
            }
            Wave::Surge { offset, peak, delay, shape } => {
                offset + peak * shape.heidler().at(t - delay)
            }
            Wave::Burst { offset, peak, delay, rate, length, period } => {
                let mut u = t - delay;
                if u < 0.0 {
                    return *offset;
                }
                if *period > 0.0 {
                    u %= period;
                }
                match burst_pulse(u, *rate, *length) {
                    Some(since) => offset + peak * FAST.at(since),
                    None => *offset,
                }
            }
        }
    }

    pub fn breaks(&self, stop: f64) -> Vec<f64> {
        match self {
            Wave::Dc | Wave::Sine { .. } => Vec::new(),
            Wave::Pulse { delay, rise, fall, width, period, .. } => {
                let mut out = Vec::new();
                let mut start = *delay;
                loop {
                    for e in [0.0, *rise, rise + width, rise + width + fall] {
                        if start + e <= stop {
                            out.push(start + e);
                        }
                    }
                    if *period <= 0.0 || start + period > stop {
                        break;
                    }
                    start += period;
                }
                out
            }
            Wave::Steps { points } => points.iter().map(|p| p.0).filter(|&t| t <= stop).collect(),
            Wave::Surge { delay, shape, .. } => {
                let h = shape.heidler();
                let mut out = Vec::new();
                marks(*delay, h.span() * 0.05, 10.0 * h.rise.max(h.fall), stop, &mut out);
                out.retain(|t| *t <= stop);
                out
            }
            Wave::Burst { delay, rate, length, period, .. } => {
                let mut out = Vec::new();
                let mut start = *delay;
                while start <= stop {
                    let count = (length * rate).ceil() as usize;
                    for k in 0..count {
                        let at = start + k as f64 / rate;
                        if at > stop {
                            break;
                        }
                        marks(at, FAST.span() * 0.1, 10.0 * FAST.fall, stop, &mut out);
                    }
                    if *period <= 0.0 {
                        break;
                    }
                    start += period;
                }
                out.retain(|t| *t <= stop);
                out
            }
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Drive {
    pub net: String,
    pub p: Node,
    pub n: Node,
    pub current: bool,
    pub dc: f64,
    pub ac: f64,
    pub r: f64,
    pub wave: Wave,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct Sweep {
    pub from: f64,
    pub to: f64,
    pub per_decade: usize,
}

impl Sweep {
    pub fn freqs(&self) -> Vec<f64> {
        let decades = (self.to / self.from).log10().max(0.0);
        let n = ((decades * self.per_decade as f64).ceil() as usize).max(1);
        (0..=n).map(|k| self.from * 10f64.powf(decades * k as f64 / n as f64)).collect()
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Noise {
    pub output: Node,
    pub output_name: String,
    pub input: Option<usize>,
    pub sweep: Sweep,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "kind", content = "of", rename_all = "snake_case")]
pub enum Measure {
    Volts(Node),
    Amps(usize),
    Phase(Node),
    Noise(Node),
    Peak(Node),
    Trough(Node),
    Rms(Node),
    Freq(Node),
    Thd(Node),
    Harmonics { element: usize, mains: f64, class_b: bool },
    Emission { node: Node, class_a: bool, average: bool },
}

impl Measure {
    pub fn over_window(&self) -> bool {
        matches!(
            self,
            Measure::Peak(_)
                | Measure::Trough(_)
                | Measure::Rms(_)
                | Measure::Freq(_)
                | Measure::Thd(_)
                | Measure::Harmonics { .. }
                | Measure::Emission { .. }
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
#[serde(tag = "kind", content = "value", rename_all = "snake_case")]
pub enum At {
    Op,
    Time(f64),
    Freq(f64),
    Window(f64, f64),
}

#[derive(Clone, Debug, Serialize)]
pub struct Check {
    pub label: String,
    pub measure: Measure,
    pub at: At,
    pub min: Option<f64>,
    pub max: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct CircuitSpec {
    pub schematic: String,
    pub nets: Vec<String>,
    pub elements: Vec<Elem>,
    pub drives: Vec<Drive>,
    pub ac: Option<Sweep>,
    pub tran: Option<(f64, f64)>,
    pub noise: Option<Noise>,
    pub checks: Vec<Check>,
    pub record: Vec<usize>,
    pub temperature: f64,
    pub netlist_hash: u64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Reading {
    pub label: String,
    pub value: f64,
    pub unit: String,
    pub passed: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max: Option<f64>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Curve {
    pub name: String,
    pub values: Vec<f64>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CircuitResult {
    pub name: String,
    pub spec_hash: u64,
    pub op: BTreeMap<String, f64>,
    pub currents: BTreeMap<String, f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub freqs: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub gain_db: Vec<Curve>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub phase_deg: Vec<Curve>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub times: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub waves: Vec<Curve>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub spectrum_freqs: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub spectrum: Vec<Curve>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub noise_freqs: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub noise: Vec<Curve>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub noise_rms: Option<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub harmonics: Vec<Curve>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub emission_freqs: Vec<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub emission: Vec<Curve>,
    pub checks: Vec<Reading>,
    pub failures: Vec<String>,
    pub passed: bool,
    pub seconds: f64,
}

pub fn fmt_si(v: f64, unit: &str) -> String {
    let a = v.abs();
    let (m, p) = if a == 0.0 {
        (1.0, "")
    } else if a >= 1e9 {
        (1e-9, "G")
    } else if a >= 1e6 {
        (1e-6, "M")
    } else if a >= 1e3 {
        (1e-3, "k")
    } else if a >= 1.0 {
        (1.0, "")
    } else if a >= 1e-3 {
        (1e3, "m")
    } else if a >= 1e-6 {
        (1e6, "u")
    } else if a >= 1e-9 {
        (1e9, "n")
    } else if a >= 1e-12 {
        (1e12, "p")
    } else {
        (1e15, "f")
    };
    let x = v * m;
    let s = format!("{x:.4}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    format!("{s}{p}{unit}")
}

pub fn class_a(order: usize) -> Option<f64> {
    let n = order as f64;
    match order {
        2 => Some(1.08),
        3 => Some(2.30),
        4 => Some(0.43),
        5 => Some(1.14),
        6 => Some(0.30),
        7 => Some(0.77),
        9 => Some(0.40),
        11 => Some(0.33),
        13 => Some(0.21),
        8..=40 if order.is_multiple_of(2) => Some(0.23 * 8.0 / n),
        15..=39 => Some(0.15 * 15.0 / n),
        _ => None,
    }
}

pub fn cispr11(f: f64, class_a: bool, average: bool) -> Option<f64> {
    if !(150e3..=30e6).contains(&f) {
        return None;
    }
    let low =
        |from: f64, to: f64| from + (to - from) * (f / 150e3).log10() / (0.5f64 / 0.15).log10();
    Some(match (class_a, average) {
        (true, false) if f < 500e3 => 79.0,
        (true, false) => 73.0,
        (true, true) if f < 500e3 => 66.0,
        (true, true) => 60.0,
        (false, false) if f < 500e3 => low(66.0, 56.0),
        (false, false) if f <= 5e6 => 56.0,
        (false, false) => 60.0,
        (false, true) if f < 500e3 => low(56.0, 46.0),
        (false, true) if f <= 5e6 => 46.0,
        (false, true) => 50.0,
    })
}
