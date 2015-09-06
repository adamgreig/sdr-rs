extern crate num;
use num::Complex;

pub type IQ<T> = Complex<T>;
pub type Real<T> = T;

pub mod downconverter;
pub mod cic;
pub mod fir;
