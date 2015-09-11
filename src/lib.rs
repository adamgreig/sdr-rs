extern crate num;
use num::Complex;

pub type IQ<T> = Complex<T>;
pub type Real<T> = T;

pub mod downconverter;
pub use downconverter::downconvert_fs_4;

pub mod cic;
pub use cic::CIC;

pub mod fir;
pub use fir::FIR;

pub mod fm;
pub use fm::fm_demod;
