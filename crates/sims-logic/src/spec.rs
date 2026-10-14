use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Level {
    L,
    H,
    X,
    Z,
}

impl Level {
    pub fn of(b: bool) -> Level {
        if b { Level::H } else { Level::L }
    }

    pub fn bit(self) -> Option<bool> {
        match self {
            Level::L => Some(false),
            Level::H => Some(true),
            _ => None,
        }
    }

    pub fn flip(self) -> Level {
        match self {
            Level::L => Level::H,
            Level::H => Level::L,
            _ => Level::X,
        }
    }

    pub fn char(self) -> char {
        match self {
            Level::L => '0',
            Level::H => '1',
            Level::X => 'x',
            Level::Z => 'z',
        }
    }

    pub fn from_char(c: char) -> Option<Level> {
        match c {
            '0' => Some(Level::L),
            '1' => Some(Level::H),
            'x' | 'X' => Some(Level::X),
            'z' | 'Z' => Some(Level::Z),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Edge {
    #[default]
    Rising,
    Falling,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OnViolation {
    #[default]
    X,
    Keep,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum GateOp {
    And,
    Or,
    Xor,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub enum Prim {
    Const(Level),
    Gate { op: GateOp, invert: bool },
    Tri,
    Dff,
    Jk,
    Sr,
    Dlatch,
    Mux2,
    Dec138,
    Counter { sync_reset: bool },
    Shift164,
    Shift595,
    Xcvr,
    Truth(Vec<(Vec<Option<bool>>, Vec<Level>)>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Clocked {
    pub clock: usize,
    pub falling: bool,
    pub data: &'static [usize],
    pub resets: &'static [usize],
    pub lead: Option<usize>,
    pub state: std::ops::Range<usize>,
}

impl Prim {
    pub fn clocked(&self) -> Vec<Clocked> {
        let on = |clock, data, resets, state| Clocked {
            clock,
            falling: false,
            data,
            resets,
            lead: None,
            state,
        };
        match self {
            Prim::Dff => vec![on(1, &[0], &[2, 3], 0..2)],
            Prim::Jk => vec![on(2, &[0, 1], &[3, 4], 0..2)],
            Prim::Dlatch => vec![Clocked { falling: true, ..on(1, &[0], &[2], 0..2) }],
            Prim::Counter { sync_reset: false } => vec![on(1, &[2, 3, 4, 5, 6, 7, 8], &[0], 0..4)],
            Prim::Counter { sync_reset: true } => vec![on(1, &[0, 2, 3, 4, 5, 6, 7, 8], &[], 0..4)],
            Prim::Shift164 => vec![on(2, &[0, 1], &[3], 0..8)],
            Prim::Shift595 => {
                vec![on(1, &[0], &[3], 0..8), Clocked { lead: Some(1), ..on(2, &[], &[], 8..16) }]
            }
            _ => Vec::new(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Input {
    Net { net: usize, invert: bool },
    Fixed(Level),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Output {
    pub net: Option<usize>,
    pub invert: bool,
    pub delay: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Cell {
    pub part: String,
    pub prim: Prim,
    pub input_names: Vec<String>,
    pub inputs: Vec<Input>,
    pub outputs: Vec<Output>,
    pub weak: bool,
    pub open_drain: bool,
    pub setup: u64,
    pub hold: u64,
    pub recovery: u64,
    pub removal: u64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Circuit {
    pub nets: Vec<String>,
    pub cells: Vec<Cell>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub enum Wave {
    Clock { period: u64, high: u64, phase: u64 },
    Steps(Vec<(u64, Level)>),
}

#[derive(Clone, Debug, Serialize)]
pub struct Stimulus {
    pub net: usize,
    pub name: String,
    pub wave: Wave,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub enum Check {
    At {
        time: u64,
        value: Vec<Option<Level>>,
    },
    Sequence {
        clock: usize,
        clock_name: String,
        edge: Edge,
        from: u64,
        values: Vec<Vec<Option<Level>>>,
    },
}

#[derive(Clone, Debug, Serialize)]
pub struct Expect {
    pub label: String,
    pub nets: Vec<usize>,
    pub names: Vec<String>,
    pub check: Check,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bus {
    pub name: String,
    pub nets: Vec<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MarkKind {
    Assertion,
    Timing,
    Contention,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Mark {
    pub time: u64,
    pub kind: MarkKind,
    pub nets: Vec<String>,
    pub text: String,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct LogicSpec {
    pub schematic: String,
    pub duration: u64,
    pub record: Vec<(String, usize)>,
    pub buses: Vec<Bus>,
    pub on_violation: OnViolation,
    pub stimuli: Vec<Stimulus>,
    pub expects: Vec<Expect>,
    pub circuit: Circuit,
    pub netlist_hash: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Trace {
    pub name: String,
    pub times: Vec<u64>,
    pub values: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BusGroup {
    pub name: String,
    pub stem: String,
    pub msb: u32,
    pub lsb: u32,
    pub bits: Vec<usize>,
}

pub fn bus_bit(name: &str) -> Option<(String, u32, String)> {
    if let Some(stem) = name.strip_suffix(']')
        && let Some((prefix, idx)) = stem.rsplit_once('[')
        && !prefix.is_empty()
        && let Ok(i) = idx.parse()
    {
        return Some((prefix.to_string(), i, String::new()));
    }
    let end = name.rfind(|c: char| c.is_ascii_digit())? + 1;
    let start = name[..end].rfind(|c: char| !c.is_ascii_digit())? + 1;
    let suffix = &name[end..];
    if !["", "_N", "_n", "#"].contains(&suffix) {
        return None;
    }
    let i = name[start..end].parse().ok()?;
    Some((name[..start].to_string(), i, suffix.to_string()))
}

type AutoBus = ((String, String), Vec<(u32, usize)>);

pub fn bus_groups(traces: &[Trace], buses: &[Bus]) -> Vec<BusGroup> {
    let find = |n: &str| traces.iter().position(|t| t.name == n);
    let mut groups: Vec<BusGroup> = Vec::new();
    let mut used = vec![false; traces.len()];
    for b in buses {
        let bits: Option<Vec<usize>> = b.nets.iter().map(|n| find(n)).collect();
        if let Some(bits) = bits
            && !bits.is_empty()
        {
            bits.iter().for_each(|k| used[*k] = true);
            groups.push(BusGroup {
                name: b.name.clone(),
                stem: b.name.clone(),
                msb: bits.len() as u32 - 1,
                lsb: 0,
                bits,
            });
        }
    }
    let mut auto: Vec<AutoBus> = Vec::new();
    for (k, t) in traces.iter().enumerate() {
        if used[k] {
            continue;
        }
        let Some((prefix, i, suffix)) = bus_bit(&t.name) else { continue };
        let key = (prefix, suffix);
        match auto.iter_mut().find(|g| g.0 == key) {
            Some(g) => g.1.push((i, k)),
            None => auto.push((key, vec![(i, k)])),
        }
    }
    for ((prefix, suffix), mut bits) in auto {
        bits.sort_by_key(|b| std::cmp::Reverse(b.0));
        let distinct = bits.windows(2).all(|w| w[0].0 != w[1].0);
        if bits.len() < 2 || !distinct {
            continue;
        }
        let (hi, lo) = (bits[0].0, bits[bits.len() - 1].0);
        groups.push(BusGroup {
            name: format!("{prefix}[{hi}:{lo}]{suffix}"),
            stem: format!("{prefix}{suffix}"),
            msb: hi,
            lsb: lo,
            bits: bits.into_iter().map(|b| b.1).collect(),
        });
    }
    groups
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LogicResult {
    pub name: String,
    pub kind: String,
    pub spec_hash: u64,
    pub duration_ps: u64,
    pub end_ps: u64,
    pub traces: Vec<Trace>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub buses: Vec<Bus>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub marks: Vec<Mark>,
    pub passed: usize,
    pub failures: Vec<String>,
    pub readings: Vec<Reading>,
    pub events: u64,
    pub seconds: f64,
    pub vcd: String,
}

pub fn fmt_time(t: u64) -> String {
    let (div, unit) = match t {
        1_000_000_000.. => (1e9, "ms"),
        1_000_000.. => (1e6, "us"),
        1_000.. => (1e3, "ns"),
        0 => (1.0, "ns"),
        _ => (1.0, "ps"),
    };
    let s = format!("{:.3}", t as f64 / div);
    let s = s.trim_end_matches('0').trim_end_matches('.');
    format!("{s}{unit}")
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Reading {
    pub label: String,
    pub value: f64,
    pub unit: String,
    #[serde(default)]
    pub detail: String,
}
