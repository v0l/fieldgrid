mod cpu;
mod ntff;
mod run;
mod sim;

pub use cpu::run_cpu;
pub use ntff::{C, Far, Patches, patches};
#[cfg(not(target_arch = "wasm32"))]
pub use run::run_blocking;
pub use run::{Extras, Pulse, Record, run};
pub use sim::*;
