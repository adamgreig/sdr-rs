use num::{Complex};

pub trait SampleType: Copy {
    fn power(&self) -> f32;
}

macro_rules! impl_scalar_sampletype {
    ($t:ty, $max:expr) => {
        impl SampleType for $t {
            #[inline]
            fn power(&self) -> f32 {
                (*self as f32 * *self as f32) / ($max * $max)
            }
        }
    }
}

macro_rules! impl_complex_sampletype {
    ($t:ty, $max:expr) => {
        impl SampleType for Complex<$t> {
            #[inline]
            fn power(&self) -> f32 {
                ((*self).im as f32 * (*self).im as f32 +
                 (*self).re as f32 * (*self).re as f32)
                / ($max * $max)
            }
        }
    }
}

macro_rules! impl_sampletype {
    ($t:ty, $max:expr) => {
        impl_scalar_sampletype!($t, $max);
        impl_complex_sampletype!($t, $max);
    }
}

impl_sampletype!(i8,  (1u8 <<7)  as f32);
impl_sampletype!(i16, (1u16<<15) as f32);
impl_sampletype!(i32, (1u32<<31) as f32);
impl_sampletype!(i64, (1u64<<63) as f32);
impl_sampletype!(f32, 1.0);
impl_sampletype!(f64, 1.0);

/// Compute the required gain in dB required to bring the input block `x`
/// up to the reference level `r`, multiplying by sub-unity constant `a` to
/// prevent large shifts.
///
/// This is a fairly simple and slow implementation which is not suitable for
/// high data rates.
pub fn agc<T: SampleType>(x: &[T], a: f32, r: f32) -> f32 {
    let avg = x.iter().fold(0.0_f32, |s, v| s + v.power()) / x.len() as f32;
    let err = 10.0 * r.log10() - 10.0 * avg.log10();
    a * err
}

#[cfg(test)]
mod tests {
    use super::agc;
    use num::Complex;

    #[test]
    fn test_agc_f32() {
        let x: Vec<f32> = vec!{0.1, 0.1, 0.1, 0.1};
        let gain = agc(&x, 1.0, 1.0);
        assert_eq!(gain, 20.0);

        let rt2 = 2.0_f32.sqrt();
        let x: Vec<f32> = vec!{0.1 * rt2, 0.0, -0.1 * rt2, 0.0};
        let gain = agc(&x, 1.0, 1.0);
        assert_eq!(gain, 20.0);

        let x: Vec<f32> = vec!{0.5, 0.5, 0.5, 0.5};
        let gain = agc(&x, 1.0, 1.0);
        assert_eq!(gain, 6.0206003);

        let x: Vec<f32> = vec!{0.1, 0.1, 0.1, 0.1};
        let gain = agc(&x, 1.0, 0.5);
        assert_eq!(gain, 20.0 - 3.0103);
    }

    #[test]
    fn test_agc_i16() {
        let x: Vec<i16> = vec!{3277, 3277, 3277, 3276};
        let gain = agc(&x, 1.0, 1.0);
        assert_eq!(gain, 20.000134);
    }

    #[test]
    fn test_agc_complex() {
        let rt2 = 2.0_f32.sqrt();
        let x: Vec<Complex<f32>> = vec!{
            Complex::new(0.1, 0.0), Complex::new(0.0, 0.1),
            Complex::new(-0.1, 0.0), Complex::new(0.0, -0.1),
        };
        let gain = agc(&x, 1.0, 1.0);
        assert_eq!(gain, 20.0);

        let x: Vec<Complex<f32>> = vec!{
            Complex::new(0.1/rt2, 0.1/rt2),
            Complex::new(0.1/rt2, 0.1/rt2),
            Complex::new(0.1/rt2, 0.1/rt2),
            Complex::new(0.1/rt2, 0.1/rt2),
        };
        let gain = agc(&x, 1.0, 1.0);
        assert_eq!(gain, 20.0);
    }
}
