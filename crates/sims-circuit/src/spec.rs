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
