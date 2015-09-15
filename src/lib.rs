extern crate num;
use num::Complex;

pub type IQ<T> = Complex<T>;
pub type Real<T> = T;

#[macro_export]
macro_rules! chain_blocks {
    ($x:ident, $($f:ident),*) => {{
        let x = $x;
        $(
            let x = $f.process(&x);
        )*
        x
    }}
}

pub mod downconverter;
pub use downconverter::RealFs4Downconverter;

pub mod cic;
pub use cic::CIC;

pub mod fir;
pub use fir::FIR;

pub mod fm;
pub use fm::FMDemod;

pub mod agc;
pub use agc::agc;
