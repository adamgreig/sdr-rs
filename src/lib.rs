extern crate num;
use num::Complex;

pub type IQ<T> = Complex<T>;
pub type Real<T> = T;

pub mod downconverter;
pub use downconverter::Downconverter;

pub mod cic;
pub use cic::CIC;

pub mod fir;
pub use fir::FIR;

pub mod prelude {
    pub use super::{IQ, Real};
    pub use super::{Downconverter, CIC, FIR};
}
