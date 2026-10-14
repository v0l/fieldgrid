#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Air {
    pub celsius: f64,
    pub humidity: f64,
    pub kilopascal: f64,
}

impl Default for Air {
    fn default() -> Air {
        Air { celsius: 20.0, humidity: 50.0, kilopascal: 101.325 }
    }
}

const REFERENCE_KELVIN: f64 = 293.15;
const REFERENCE_KILOPASCAL: f64 = 101.325;
const TRIPLE_KELVIN: f64 = 273.16;

impl Air {
    pub fn kelvin(&self) -> f64 {
        self.celsius + 273.15
    }

    pub fn speed(&self) -> f64 {
        343.2 * (self.kelvin() / REFERENCE_KELVIN).sqrt()
    }

    pub fn decibels_per_metre(&self, hz: f64) -> f64 {
        let t = self.kelvin();
        let pa = self.kilopascal / REFERENCE_KILOPASCAL;
        let saturation =
            REFERENCE_KILOPASCAL * 10f64.powf(-6.8346 * (TRIPLE_KELVIN / t).powf(1.261) + 4.6151);
        let h = self.humidity * saturation / self.kilopascal;
        let tr = t / REFERENCE_KELVIN;
        let oxygen = pa * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
        let nitrogen =
            pa * tr.powf(-0.5) * (9.0 + 280.0 * h * (-4.170 * (tr.powf(-1.0 / 3.0) - 1.0)).exp());
        8.686
            * hz
            * hz
            * (1.84e-11 / pa * tr.sqrt()
                + tr.powf(-2.5)
                    * (0.01275 * (-2239.1 / t).exp() / (oxygen + hz * hz / oxygen)
                        + 0.1068 * (-3352.0 / t).exp() / (nitrogen + hz * hz / nitrogen)))
    }

    pub fn nepers_per_metre(&self, hz: f64) -> f64 {
        self.decibels_per_metre(hz) * std::f64::consts::LN_10 / 10.0
    }
}
