use num::{Num, NumCast, Zero, Complex};
use std::f64::consts::PI;

/// Implement the underlying operations and types required by the FIR.
///
/// In general the AccType should be larger than the SampleType but compatible
/// (i.e. you can add a SampleType to an AccType), while the TapType should
/// always be scalar and is usually the same size as the AccType.
///
/// gain() specifies the gain of the filter, i.e. the sum of all the taps,
/// divided by the interpolation level L.
/// input() copies an item from the input into a given delay line element.
/// accumulate() updates an accumulator with a delay line element * a tap.
/// output() writes to the output array from an accumulator.
pub trait SampleType: Copy {
    type AccType: Zero + Copy;
    type TapType: Num + NumCast + Zero + Copy;
    fn gain(interpolate: usize) -> Self::TapType;
    unsafe fn input(x: *const Self, delay: *mut Self::AccType);
    unsafe fn accumulate(acc: &mut Self::AccType, delay: *const Self::AccType,
                         tap: *const Self::TapType);
    unsafe fn output(acc: Self::AccType, out: *mut Self, gain: Self::TapType);
}

/// Implement SampleType for a scalar type such as i16 or f32.
/// $t is the sample type, $tt the acc/tap type and $tapsum the filter gain.
macro_rules! impl_scalar_sampletype {
    ($t:ty, $tt:ty, $tapsum:expr) => {
        impl SampleType for $t {
            type AccType = $tt;
            type TapType = $tt;
            
            #[inline]
            fn gain(interpolate: usize) -> $tt { $tapsum / interpolate as $tt }
            
            #[inline]
            unsafe fn input(x: *const $t, delay: *mut $tt) {
                *delay = *x as $tt;
            }

            #[inline]
            unsafe fn accumulate(acc: &mut $tt, delay: *const $tt,
                                 tap: *const $tt) {
                *acc += *delay * *tap;
            }

            #[inline]
            unsafe fn output(acc: $tt, out: *mut $t, gain: $tt) {
                *out = (acc / gain) as $t
            }
        }
    }
}

/// Implement SampleType for a Complex type such as Complex<i16>.
/// $t is the sample type, $tt the acc/tap type and $tapsum the filter gain.
macro_rules! impl_complex_sampletype {
    ($t:ty, $tt:ty, $tapsum:expr) => {
        impl SampleType for Complex<$t> {
            type AccType = Complex<$tt>;
            type TapType = $tt;

            #[inline]
            fn gain(interpolate: usize) -> $tt { $tapsum / interpolate as $tt }

            #[inline]
            unsafe fn input(x: *const Complex<$t>, delay: *mut Complex<$tt>) {
                (*delay).re = (*x).re as $tt;
                (*delay).im = (*x).im as $tt;
            }

            #[inline]
            unsafe fn accumulate(acc: &mut Complex<$tt>,
                                 delay: *const Complex<$tt>, tap: *const $tt) {
                (*acc).re += (*delay).re * *tap;
                (*acc).im += (*delay).im * *tap;
            }

            #[inline]
            unsafe fn output(acc: Complex<$tt>, out: *mut Complex<$t>,
                             gain: $tt) {
                (*out).re = (acc.re / gain) as $t;
                (*out).im = (acc.im / gain) as $t;
            }
        }
    }
}

/// Implement Scalar and Complex SampleTypes for the same underlying types.
macro_rules! impl_sampletype {
    ($t:ty, $tt:ty, $tapsum:expr) => {
        impl_scalar_sampletype!($t, $tt, $tapsum);
        impl_complex_sampletype!($t, $tt, $tapsum);
    }
}

impl_sampletype!(i8,  i16, 1<<7 );
impl_sampletype!(i16, i32, 1<<15);
impl_sampletype!(i32, i64, 1<<31);
impl_sampletype!(f32, f64, 1.0);

/// FIR filter.
pub struct FIR<T: SampleType> {
    taps: Vec<T::TapType>,
    tap_idx: isize,
    delay: Vec<T::AccType>,
    delay_idx: isize,
    decimate: usize,
    interpolate: usize,
}

impl <T: SampleType> FIR<T> {
    /// Create a new FIR with the given taps and decimation.
    ///
    /// Taps should sum to T::gain() or close to it.
    ///
    /// Set decimate=1 for no decimation, decimate=2 for /2, etc.
    /// Set interpolate=1 for no interpolation, interplate=2 for *2, etc.
    /// Note that if the number of taps is not a multiple of the interpolation
    /// ratio then they will be zero-padded at the end until they are.
    ///
    /// Implements a polyphase FIR to do efficient decimation and interpolation
    /// (identical to a standard FIR when interpolate=1).
    pub fn new(taps: &Vec<T::TapType>, decimate: usize, interpolate: usize)
        -> FIR<T>
    {
        assert!(taps.len() > 0);

        // Copy the taps and zero-pad to get a multiple of interpolation ratio
        let mut taps: Vec<T::TapType> = taps.to_owned();
        if taps.len() % interpolate != 0 {
            for _ in 0..(interpolate - (taps.len() % interpolate)) {
                taps.push(T::TapType::zero());
            }
        }

        // Set up and fill the delay line with zeros
        let n_delay = taps.len() / interpolate;
        let mut delay: Vec<T::AccType> = Vec::with_capacity(n_delay);
        for _ in 0..n_delay {
            delay.push(T::AccType::zero());
        }

        FIR { taps: taps, tap_idx: interpolate as isize - 1isize,
              delay: delay, delay_idx: 0isize,
              decimate: decimate, interpolate: interpolate }
    }

    /// Create a new FIR from a number of taps, desired frequency response
    /// (a 512-long vector from 0 to Nyquist freq) and decimation.
    ///
    /// Set decimate=1 for no decimation, decimate=2 for /2, etc.
    /// Set interpolate=1 for no interpolation, interplate=2 for *2, etc.
    pub fn from_gains(n_taps: usize, gains: &Vec<f64>,
                      decimate: usize, interpolate: usize)
        -> FIR<T>
    {
        let taps = firwin2(n_taps, gains);
        let taps = quantise_taps::<T::TapType>(&taps, T::gain(1));
        FIR::new(&taps, decimate, interpolate)
    }

    /// Create a new FIR that compensates for a CIC filter specified by
    /// q and r, optionally also cutting off the frequency response after
    /// Fs/(2*decimate) and having the FIR decimate by that factor.
    ///
    /// Set decimate=1 for no decimation, decimate=2 for /2, etc.
    ///
    /// TODO: Does not quite match results obtained in Python, with slightly
    /// worse simulated performance. Investigate.
    pub fn cic_compensator(n_taps: usize, q: usize, r: usize, decimate: usize)
        -> FIR<T>
    {
        let q = q as i32;
        let r = r as f64;
        let f: Vec<f64> = (0..512).map(|x|
            x as f64 / (16.0_f64 * 512.0_f64)).collect();
        let mut h: Vec<f64> = f.iter().map(|f|
            ((PI * f * r).sin() / (r * (PI * f).sin())).abs().powi(q).recip())
            .collect();
        h[0] = 1.0_f64;
        h[511] = 0.0_f64;
        if decimate > 1 {
            for n in (512/decimate - 51/decimate)..512 {
                h[n] = 0.0_f64;
            }
        }
        FIR::from_gains(n_taps, &h, decimate, 1) 
    }

    /// Create a new FIR that implements a specified raised cosine filter,
    /// with rolloff factor b, reciprocal symbol rate t in normalised time
    /// (i.e. samples per symbol), and half-width w (in bit periods, i.e. how
    /// many periods of t to extend either side of centre).
    pub fn from_rcf(b: f64, t: usize, w: usize) -> FIR<T> {
        assert!(b >= 0.0_f64 && b <= 1.0_f64);
        let t = t as isize;
        let w = w as isize;
        let n_taps = 2 * t * w + 1;
        let mut taps: Vec<f64> = Vec::with_capacity(n_taps as usize);
        for i in (-t*w)..(t*w + 1) {
            let x = i as f64 / t as f64;
            let mut y: f64;

            if i == 0 {
                y = 1.0;
            } else if i % t == 0 {
                y = 0.0;
            } else {
                y = (PI * x).sin() / (PI * x);
                y *= (PI * b * x).cos();
                y /= 1.0 - (4.0 * b*b * x*x);
            }
            taps.push(y);
        }
        let taps = quantise_taps::<T::TapType>(&taps, T::gain(1));
        FIR::new(&taps, 1, 1)
    }

    /// Create a new FIR resampler, with frequency response suitable for the
    /// resampling ratio. n_taps should ideally be a multiple of interpolate
    /// (but will be zero-padded at the end if not).
    pub fn resampler(n_taps: usize, decimate: usize, interpolate: usize)
        -> FIR<T>
    {
        let fc = 512 / ::std::cmp::max(decimate, interpolate);
        let mut gains: Vec<f64> = Vec::with_capacity(512);
        for _ in 0..fc {
            gains.push(1.0);
        }
        for _ in fc..512 {
            gains.push(0.0);
        }
        FIR::from_gains(n_taps, &gains, decimate, interpolate)
    }

    /// Return a reference to the filter's taps.
    pub fn taps(&self) -> &Vec<T::TapType> {
        &self.taps
    }

    /// Process a block of data x, outputting the filtered and possibly
    /// decimated data.
    pub fn process(&mut self, x: &Vec<T>) -> Vec<T> {
        // Check we were initialised correctly and
        // ensure invariances required for unsafe code.
        assert!(self.taps.len() % self.interpolate == 0);
        assert!(self.delay.len() == self.taps.len() / self.interpolate);
        assert!(self.delay_idx < self.delay.len() as isize);
        assert!(self.tap_idx < (self.taps.len() / self.interpolate) as isize);
        assert_eq!(x.len() % self.decimate, 0);

        // Allocate output
        let mut y: Vec<T> = Vec::with_capacity(
            (x.len() * self.interpolate) / self.decimate);
        unsafe { y.set_len((x.len() * self.interpolate) / self.decimate) };

        // Grab pointers to various things
        let mut delay_idx = self.delay_idx;
        let mut tap_idx = self.tap_idx;
        let delay_len = self.delay.len() as isize;
        let delay_p = &mut self.delay[0] as *mut T::AccType;
        let out_p = &mut y[0] as *mut T;
        let mut in_p = &x[0] as *const T;
        let ylen = y.len() as isize;
        let decimate = self.decimate as isize;
        let interpolate = self.interpolate as isize;
        let tap0 = &self.taps[0] as *const T::TapType;
        let gain: T::TapType = T::gain(self.interpolate);

        // Process each actually generated output sample
        for k in 0..ylen {

            // For every high-rate clock tick, advance the polyphase
            // coefficient commutators by one, and when they wrap around,
            // insert a new input into the delay line and advance that by one.
            // Repeat until a sample we're not going to skip.
            for _ in 0..decimate {
                // Advance coefficient commutators.
                // Note that the initialised value for tap_idx is
                // interpolate - 1, so that on the first run we'll reset it
                // to zero and add the first sample to the delay line.
                tap_idx += 1;
                if tap_idx == interpolate {
                    tap_idx = 0;

                    // Insert input sample
                    unsafe {
                        T::input(in_p, delay_p.offset(delay_idx));
                        in_p = in_p.offset(1);
                    }
                    delay_idx -= 1;
                    if delay_idx == -1 {
                        delay_idx = delay_len - 1;
                    }
                }
            }

            // Compute the multiply-accumulate only for output samples.
            let mut acc: T::AccType = T::AccType::zero();
            let mut tap_p = unsafe { tap0.offset(tap_idx) };

            // First compute from the index to the end of the buffer
            for idx in (delay_idx + 1)..delay_len {
                unsafe {
                    T::accumulate(&mut acc, delay_p.offset(idx), tap_p);
                    tap_p = tap_p.offset(interpolate);
                }
            }

            // Then compute from the start of the buffer to the index
            for idx in 0..(delay_idx + 1) {
                unsafe {
                    T::accumulate(&mut acc, delay_p.offset(idx), tap_p);
                    tap_p = tap_p.offset(interpolate);
                }
            }

            // Save the result, accounting for filter gain
            unsafe { T::output(acc, out_p.offset(k), gain); }
        }

        // Update index for next time
        self.delay_idx = delay_idx;
        self.tap_idx = tap_idx;

        y
    }
}


/// Design an FIR filter using the window method.
///
/// n_taps: the number of taps to return
/// gains: desired frequency response evalulated on a 512-point grid between
///        zero and the Nyquist frequency
pub fn firwin2(n_taps: usize, gains: &Vec<f64>) -> Vec<f64> {
    assert!(n_taps > 0);
    assert_eq!(gains.len(), 512);
    assert_eq!(gains[511], 0.0f64);

    // Gather complex gains
    let mut gainsc: Vec<Complex<f64>> =
        gains.iter().map(|g| Complex::new(*g, 0.0)).collect();

    // Phase shift gains so the resulting impulse response is nicely aligned
    // to just grab the first n_taps coefficients.
    for (idx, gain) in gainsc.iter_mut().enumerate() {
        let (r, mut p) = gain.to_polar();
        p += -(n_taps as f64 - 1.0)/2.0 * PI * (idx as f64)/512.0;
        *gain = Complex::from_polar(&r, &p);
    }

    // Inverse the frequency response to get the impulse response (aka taps)
    let taps: Vec<f64> = irdft(&gainsc).into_iter().take(n_taps).collect();

    // Multiply by window
    let w = hamming(n_taps);
    taps.iter().zip(w.iter()).map(|(t, w)| t * w).collect()
}

/// Quantise FIR taps to i16 taps that sum to `total`
pub fn quantise_taps<T: Num + NumCast>(taps: &Vec<f64>, total: T) -> Vec<T> {
    let sum: f64 = taps.iter().fold(0.0_f64, |acc, &x| acc + x);
    let total: f64 = <f64 as NumCast>::from(total).unwrap();
    taps.iter().map(|t| <T as NumCast>::from(t * total / sum).unwrap()).collect()
}

/// Compute the Hamming window over n samples
fn hamming(n: usize) -> Vec<f64> {
    (0..n).map(|x| 0.54 - 0.46*(2.0 * PI * x as f64 / (n as f64 - 1.0)).cos())
          .collect()
}

/// Compute the DFT of x.
#[allow(non_snake_case)]
fn dft(x: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let N = x.len();
    let mut out: Vec<Complex<f64>> = Vec::with_capacity(N);
    for k in 0..N {
        out.push(Complex::new(0.0, 0.0));
        for n in 0..N {
            let f = 2.0 * PI * (k as f64) * (n as f64) / (N as f64);
            out[k] = out[k] + x[n] * Complex::new(f.cos(), -f.sin());
        }
    }
    out
}

/// Compute the IDFT of x.
///
/// IDFT(x) = conj(DFT(conj(x))) / N
#[allow(non_snake_case)]
fn idft(x: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let N = x.len();
    let x = x.iter().map(|x| x.conj()).collect();
    let y = dft(&x);
    y.iter().map(|y| y.conj().unscale(N as f64)).collect()
}

/// Compute the IRDFT of x.
#[allow(non_snake_case)]
fn irdft(x: &Vec<Complex<f64>>) -> Vec<f64> {
    let No2p1 = x.len();
    let No2 = No2p1 - 1;
    let mut xc = x.to_owned();
    for n in 1..No2 {
        xc.push(x[No2 - n].conj());
    }
    idft(&xc).iter().map(|x| x.re).collect()
}

#[cfg(test)]
mod tests {
    use num::Complex;
    use super::{FIR, firwin2, quantise_taps, hamming, dft, idft, irdft};

    #[test]
    fn test_fir_impulse() {
        let taps: Vec<i32> = vec!{8192, 16384, 8192};
        let mut fir = FIR::new(&taps, 1, 1);
        let x: Vec<i16> = vec!{4, 0, 0, 0, 0, 0, 0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{1, 2, 1, 0, 0, 0, 0});
    }

    #[test]
    fn test_fir_impulse_f() {
        let taps: Vec<f64> = vec!{0.25, 0.5, 0.25};
        let mut fir = FIR::new(&taps, 1, 1);
        let x: Vec<f32> = vec!{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{0.25, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0});
    }

    #[test]
    fn test_fir_impulse_c() {
        let taps: Vec<i32> = vec!{8192, 16384, 8192};
        let mut fir = FIR::new(&taps, 1, 1);
        let x: Vec<Complex<i16>> = vec!{
            Complex{re:4, im:-8},
            Complex{re:0, im: 0},
            Complex{re:0, im: 0},
            Complex{re:0, im: 0},
            Complex{re:0, im: 0},
        };
        let y = fir.process(&x);
        assert_eq!(y, vec!{
            Complex{re:1, im:-2},
            Complex{re:2, im:-4},
            Complex{re:1, im:-2},
            Complex{re:0, im: 0},
            Complex{re:0, im: 0}});
    }

    #[test]
    fn test_fir_impulse_cf() {
        let taps: Vec<f64> = vec!{0.25, 0.5, 0.25};
        let mut fir = FIR::new(&taps, 1, 1);
        let x: Vec<Complex<f32>> = vec!{
            Complex{re:1.0, im:-2.0},
            Complex{re:0.0, im: 0.0},
            Complex{re:0.0, im: 0.0},
            Complex{re:0.0, im: 0.0},
            Complex{re:0.0, im: 0.0},
        };
        let y = fir.process(&x);
        assert_eq!(y, vec!{
            Complex{re:0.25, im:-0.5},
            Complex{re:0.50, im:-1.0},
            Complex{re:0.25, im:-0.5},
            Complex{re:0.00, im: 0.0},
            Complex{re:0.00, im: 0.0}});
    }

    #[test]
    fn test_fir_decimate() {
        let taps: Vec<i32> = vec!{8192, 8192, 8192, 8192};
        let mut fir = FIR::new(&taps, 2, 1);
        let x: Vec<i16> = vec!{4, 4, 4, 4, 8, 8, 8, 8};
        let y = fir.process(&x);
        assert_eq!(y, vec!{2, 4, 6, 8});
    }

    #[test]
    fn test_fir_interpolate() {
        let taps: Vec<i32> = vec!{
            -55, 0, 96, 0, -220, 0, 461, 0, -877, 0, 1608, 0, -3176, 0, 10342,
            16410, 10342, 0, -3176, 0, 1608, 0, -877, 0, 461, 0, -220, 0, 96,
            0, -55};
        let mut fir = FIR::new(&taps, 1, 2);
        let x: Vec<i16> = vec!{
            10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
            160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280,
            290, 300};
        let y = fir.process(&x);
        assert_eq!(y, vec!{
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 3, 10, 15, 20, 24, 30,
            34, 40, 44, 50, 54, 60, 64, 70, 74, 80, 84, 90, 94, 100, 104, 110,
            114, 120, 124, 130, 134, 140, 144, 150, 154, 160, 164, 170, 174,
            180, 184, 190, 194, 200, 204, 210, 214, 220, 224, 230});
    }

    #[test]
    fn test_fir_continuity() {
        let taps: Vec<i32> = vec!{8192, 16384, 8192};
        let mut fir = FIR::new(&taps, 1, 1);
        let x: Vec<i16> = vec!{4, 0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{1, 2});
        let x: Vec<i16> = vec!{0, 0, 0, 0, 0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{1, 0, 0, 0, 0});
    }

    #[test]
    fn test_fir_compensate() {
        let fir = FIR::<i16>::cic_compensator(63, 5, 8, 2);
        let taps = fir.taps();
        assert_eq!(*taps, vec!{
            -4, -43, -9, 53, 33, -66, -74, 72, 138, -54, -221, -3, 313, 118,
            -389, -303, 419, 560, -359, -880, 161, 1240, 241, -1602, -949,
            1899, 2182, -1975, -4636, 998, 11948, 17090, 10707, -142, -4690,
            -1379, 2428, 1591, -1214, -1464, 473, 1206, -20, -906, -230, 615,
            338, -366, -347, 176, 297, -48, -222, -23, 147, 53, -85, -56, 42,
            49, -17, -42, 2});
    }

    #[test]
    fn test_fir_rcf() {
        let fir = FIR::<i16>::from_rcf(0.5, 5, 3);
        let taps = fir.taps();
        assert_eq!(*taps, vec!{0i32, 19, 78, 140, 138, 0, -290, -644, -870,
            -719, 0, 1319, 3045, 4790, 6090, 6571, 6090, 4790, 3045, 1319, 0,
            -719, -870, -644, -290, 0, 138, 140, 78, 19, 0});
    }

    #[test]
    fn test_fir_resampler() {
        let mut fir = FIR::<i16>::resampler(20, 2, 3);
        let x = vec!{5, 10, 15, 20, 25, 30, 35, 40, 40, 40, 40, 40};
        let y = fir.process(&x);
        assert_eq!(y, vec!{
            0, 0, 0, 0, 4, 7, 10, 14, 17, 20, 24, 27, 30, 34, 38, 40, 40, 39});
    }

    #[test]
    fn test_firwin2() {
        let mut gains: Vec<f64> = Vec::with_capacity(512);
        for i in 0..512 {
            if i < 128 {
                gains.push(1.0_f64);
            } else {
                gains.push(0.0_f64);
            }
        }
        let taps = firwin2(16, &gains);
        assert_eq!(taps, vec!{
            -0.0013739594312205646, -0.005478372904749667,
            -0.012343994610893342, -0.010323351794476321, 0.021387446454536923,
            0.09167900348875416, 0.17954532405889292, 0.24108715546226311,
            0.24035744268642953, 0.17776648555261024, 0.08976231035333902,
            0.020081701014194094, -0.01085320791569337, -0.012398611604538569,
            -0.00540355168901675, -0.0012970506629601231});
    }

    #[test]
    fn test_quantise_taps() {
        let taps = vec!{0.25, 0.5, 0.25};
        assert_eq!(quantise_taps(&taps, 32768), vec!{8192, 16384, 8192});
    }

    #[test]
    fn test_hamming() {
        let h3 = hamming(3);
        let h5 = hamming(5);
        assert!((h3[0] - 0.08).abs() < 1e-8);
        assert!((h3[1] - 1.00).abs() < 1e-8);
        assert!((h3[2] - 0.08).abs() < 1e-8);
        assert!((h5[0] - 0.08).abs() < 1e-8);
        assert!((h5[1] - 0.54).abs() < 1e-8);
        assert!((h5[2] - 1.00).abs() < 1e-8);
        assert!((h5[3] - 0.54).abs() < 1e-8);
        assert!((h5[4] - 0.08).abs() < 1e-8);
    }

    #[test]
    fn test_dft_1() {
        let x = vec!{
            Complex::new(1.0, 1.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        };
        let y = dft(&x);
        assert_eq!(y, vec!{
            Complex::new(1.0, 1.0),
            Complex::new(1.0, 1.0),
            Complex::new(1.0, 1.0),
            Complex::new(1.0, 1.0),
        });
    }

    #[test]
    fn test_dft_2() {
        let x = vec!{
            Complex::new(0.0,  1.0),
            Complex::new(0.0, -1.0),
            Complex::new(0.0,  1.0),
            Complex::new(0.0, -1.0),
        };
        let y = dft(&x);
        assert!((y[0].re - 0.0).abs() < 1e-8);
        assert!((y[1].re - 0.0).abs() < 1e-8);
        assert!((y[2].re - 0.0).abs() < 1e-8);
        assert!((y[3].re - 0.0).abs() < 1e-8);
        assert!((y[0].im - 0.0).abs() < 1e-8);
        assert!((y[1].im - 0.0).abs() < 1e-8);
        assert!((y[2].im - 4.0).abs() < 1e-8);
        assert!((y[3].im - 0.0).abs() < 1e-8);
    }

    #[test]
    fn test_dft_3() {
        let x = vec!{
            Complex::new( 2.0, 1.0),
            Complex::new( 0.0, 1.0),
            Complex::new(-2.0, 1.0),
            Complex::new( 0.0, 1.0),
            Complex::new( 2.0, 1.0),
            Complex::new( 0.0, 1.0),
            Complex::new(-2.0, 1.0),
            Complex::new( 0.0, 1.0),
        };
        let y = dft(&x);
        assert!((y[0].re - 0.0).abs() < 1e-7);
        assert!((y[0].im - 8.0).abs() < 1e-7);
        assert!((y[1].re - 0.0).abs() < 1e-7);
        assert!((y[1].im - 0.0).abs() < 1e-7);
        assert!((y[2].re - 8.0).abs() < 1e-7);
        assert!((y[2].im - 0.0).abs() < 1e-7);
        assert!((y[3].re - 0.0).abs() < 1e-7);
        assert!((y[3].im - 0.0).abs() < 1e-7);
        assert!((y[4].re - 0.0).abs() < 1e-7);
        assert!((y[4].im - 0.0).abs() < 1e-7);
        assert!((y[5].re - 0.0).abs() < 1e-7);
        assert!((y[5].im - 0.0).abs() < 1e-7);
        assert!((y[6].re - 8.0).abs() < 1e-7);
        assert!((y[6].im - 0.0).abs() < 1e-7);
        assert!((y[7].re - 0.0).abs() < 1e-7);
        assert!((y[7].im - 0.0).abs() < 1e-7);
    }

    #[test]
    fn test_dft_4() {
        let x = vec!{
            Complex::new( 1.0, 0.0),
            Complex::new(-1.0, 0.0),
            Complex::new( 2.0, 0.0),
            Complex::new(-2.0, 0.0),
        };
        let y = dft(&x);
        assert!((y[0].re -  0.0).abs() < 1e-8);
        assert!((y[0].im -  0.0).abs() < 1e-8);
        assert!((y[1].re - -1.0).abs() < 1e-8);
        assert!((y[1].im - -1.0).abs() < 1e-8);
        assert!((y[2].re -  6.0).abs() < 1e-8);
        assert!((y[2].im -  0.0).abs() < 1e-8);
        assert!((y[3].re - -1.0).abs() < 1e-8);
        assert!((y[3].im -  1.0).abs() < 1e-8);
    }

    #[test]
    fn test_idft_1() {
        let x = vec!{
            Complex::new(1.0, 1.0),
            Complex::new(1.0, 1.0),
            Complex::new(1.0, 1.0),
            Complex::new(1.0, 1.0),
        };
        let y = idft(&x);
        assert!((y[0].re - 1.0).abs() < 1e-8);
        assert!((y[1].re - 0.0).abs() < 1e-8);
        assert!((y[2].re - 0.0).abs() < 1e-8);
        assert!((y[3].re - 0.0).abs() < 1e-8);
        assert!((y[0].im - 1.0).abs() < 1e-8);
        assert!((y[1].im - 0.0).abs() < 1e-8);
        assert!((y[2].im - 0.0).abs() < 1e-8);
        assert!((y[3].im - 0.0).abs() < 1e-8);
    }

    #[test]
    fn test_idft_2() {
        let x = vec!{
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 4.0),
            Complex::new(0.0, 0.0),
        };
        let y = idft(&x);
        assert!((y[0].re - 0.0).abs() < 1e-8);
        assert!((y[1].re - 0.0).abs() < 1e-8);
        assert!((y[2].re - 0.0).abs() < 1e-8);
        assert!((y[3].re - 0.0).abs() < 1e-8);
        assert!((y[0].im -  1.0).abs() < 1e-8);
        assert!((y[1].im - -1.0).abs() < 1e-8);
        assert!((y[2].im -  1.0).abs() < 1e-8);
        assert!((y[3].im - -1.0).abs() < 1e-8);
    }

    #[test]
    fn test_idft_3() {
        let x = vec!{
            Complex::new( 0.0,  0.0),
            Complex::new(-1.0, -1.0),
            Complex::new( 6.0,  0.0),
            Complex::new(-1.0,  1.0),
        };
        let y = idft(&x);
        assert!((y[0].re -  1.0).abs() < 1e-8);
        assert!((y[1].re - -1.0).abs() < 1e-8);
        assert!((y[2].re -  2.0).abs() < 1e-8);
        assert!((y[3].re - -2.0).abs() < 1e-8);
        assert!((y[0].im -  0.0).abs() < 1e-8);
        assert!((y[1].im -  0.0).abs() < 1e-8);
        assert!((y[2].im -  0.0).abs() < 1e-8);
        assert!((y[3].im -  0.0).abs() < 1e-8);
    }

    #[test]
    fn test_irdft_1() {
        let x = vec!{
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
        };
        let y = irdft(&x);
        assert!((y[0] - 1.0).abs() < 1e-8);
        assert!((y[1] - 0.0).abs() < 1e-8);
        assert!((y[2] - 0.0).abs() < 1e-8);
        assert!((y[3] - 0.0).abs() < 1e-8);
    }

    #[test]
    fn test_irdft_2() {
        let x = vec!{
            Complex::new( 0.0,  0.0),
            Complex::new(-1.0, -1.0),
            Complex::new( 6.0,  0.0),
        };
        let y = irdft(&x);
        assert!((y[0] -  1.0).abs() < 1e-8);
        assert!((y[1] - -1.0).abs() < 1e-8);
        assert!((y[2] -  2.0).abs() < 1e-8);
        assert!((y[3] - -2.0).abs() < 1e-8);
    }
}
