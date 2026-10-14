use crate::BANDS;

#[derive(Clone, Debug, PartialEq)]
pub struct Response {
    pub bin: f64,
    pub direct: f64,
    pub energy: [Vec<f64>; BANDS.len()],
}

pub const MODULATION: [f64; 14] =
    [0.63, 0.8, 1.0, 1.25, 1.6, 2.0, 2.5, 3.15, 4.0, 5.0, 6.3, 8.0, 10.0, 12.5];
const MALE_WEIGHT: [f64; 7] = [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173];
const MALE_REDUNDANCY: [f64; 6] = [0.085, 0.078, 0.065, 0.011, 0.047, 0.095];

fn fit(curve: &[f64], bin: f64, from: f64, to: f64) -> Option<f64> {
    let points: Vec<(f64, f64)> = curve
        .iter()
        .enumerate()
        .filter(|(_, v)| **v <= from && **v >= to)
        .map(|(i, v)| ((i as f64 + 0.5) * bin, *v))
        .collect();
    if points.len() < 3 {
        return None;
    }
    let n = points.len() as f64;
    let (sx, sy) = points.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0, a.1 + p.1));
    let (mx, my) = (sx / n, sy / n);
    let (sxy, sxx) = points
        .iter()
        .fold((0.0, 0.0), |a, p| (a.0 + (p.0 - mx) * (p.1 - my), a.1 + (p.0 - mx) * (p.0 - mx)));
    let slope = sxy / sxx;
    (slope < 0.0).then(|| -60.0 / slope)
}

impl Response {
    pub fn decay(&self, band: usize) -> Vec<f64> {
        let e = &self.energy[band];
        let total: f64 = e.iter().sum();
        let mut left = total;
        e.iter()
            .map(|v| {
                let level = 10.0 * (left.max(1e-300) / total).log10();
                left -= v;
                level
            })
            .collect()
    }

    pub fn reverberation(&self, band: usize, range: f64) -> Option<f64> {
        fit(&self.decay(band), self.bin, -5.0, -5.0 - range)
    }

    pub fn t20(&self, band: usize) -> Option<f64> {
        self.reverberation(band, 20.0)
    }

    pub fn t30(&self, band: usize) -> Option<f64> {
        self.reverberation(band, 30.0)
    }

    pub fn edt(&self, band: usize) -> Option<f64> {
        fit(&self.decay(band), self.bin, 0.0, -10.0)
    }

    fn split(&self, band: usize, after: f64) -> (f64, f64) {
        let cut = self.direct + after;
        let mut early = 0.0;
        let mut late = 0.0;
        for (i, v) in self.energy[band].iter().enumerate() {
            let (start, end) = (i as f64 * self.bin, (i + 1) as f64 * self.bin);
            let part = ((cut - start) / (end - start)).clamp(0.0, 1.0);
            early += v * part;
            late += v * (1.0 - part);
        }
        (early, late)
    }

    pub fn clarity(&self, band: usize, ms: f64) -> f64 {
        let (early, late) = self.split(band, ms / 1000.0);
        10.0 * (early / late).log10()
    }

    pub fn definition(&self, band: usize) -> f64 {
        let (early, late) = self.split(band, 0.05);
        early / (early + late)
    }

    pub fn modulation(&self, band: usize, hz: f64) -> f64 {
        let w = 2.0 * std::f64::consts::PI * hz;
        let (mut re, mut im, mut total) = (0.0, 0.0, 0.0);
        for (i, v) in self.energy[band].iter().enumerate() {
            let t = (i as f64 + 0.5) * self.bin;
            re += v * (w * t).cos();
            im += v * (w * t).sin();
            total += v;
        }
        (re * re + im * im).sqrt() / total
    }

    pub fn sti(&self) -> f64 {
        let mti: Vec<f64> = (0..BANDS.len())
            .map(|band| {
                MODULATION
                    .iter()
                    .map(|f| {
                        let m = self.modulation(band, *f).min(1.0 - 1e-12);
                        let snr = (10.0 * (m / (1.0 - m)).log10()).clamp(-15.0, 15.0);
                        (snr + 15.0) / 30.0
                    })
                    .sum::<f64>()
                    / MODULATION.len() as f64
            })
            .collect();
        let weighted: f64 = mti.iter().zip(MALE_WEIGHT).map(|(m, a)| m * a).sum();
        let redundant: f64 =
            mti.windows(2).zip(MALE_REDUNDANCY).map(|(m, b)| b * (m[0] * m[1]).sqrt()).sum();
        (weighted - redundant).min(1.0)
    }
}
