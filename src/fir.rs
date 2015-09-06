use num::Complex;
use std::f64::consts::PI;

/// FIR filter with i16 taps and i32 registers and optional decimation.
pub struct FIR {
    taps: Vec<i32>,
    delay: Vec<i32>,
    decimate: usize,
    delay_idx: isize
}

impl FIR {
    /// Create a new FIR with the given taps and decimation.
    ///
    /// Taps should sum to 32768 or close to it.
    ///
    /// Set decimate=1 for no decimation, decimate=2 for /2, etc.
    pub fn new(taps: &Vec<i16>, decimate: usize)
        -> FIR
    {
        assert!(taps.len() > 0);
        let taps: Vec<i32> = taps.iter().map(|t| *t as i32).collect();
        let mut delay: Vec<i32> = Vec::with_capacity(taps.len());
        for _ in 0..taps.len() {
            delay.push(0i32);
        }
        FIR { taps: taps, delay: delay,
              decimate: decimate, delay_idx: 0isize }
    }

    /// Create a new FIR from a number of taps, desired frequency response
    /// (a 512-long vector from 0 to Nyquist freq) and decimation.
    ///
    /// Set decimate=1 for no decimation, decimate=2 for /2, etc.
    pub fn from_gains(n_taps: usize, gains: &Vec<f64>, decimate: usize)
        -> FIR
    {
        FIR::new(&quantise_taps_i16(&firwin2(n_taps, gains)), decimate)
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
        -> FIR 
    {
        assert!(decimate > 0);
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
        FIR::from_gains(n_taps, &h, decimate) 
    }

    /// Return a reference to the filter's taps.
    pub fn taps(&self) -> &Vec<i32> {
        &self.taps
    }

    /// Process a block of data x, outputting the filtered and possibly
    /// decimated data.
    pub fn process(&mut self, x: &Vec<i16>) -> Vec<i16> {
        // Check we were initialised correctly and
        // ensure invariances required for unsafe code.
        assert!(self.taps.len() > 0);
        assert!(self.decimate > 0);
        assert!(self.decimate < self.delay.len());
        assert!(self.delay_idx < self.delay.len() as isize);
        assert_eq!(self.delay.len(), self.taps.len());
        assert_eq!(x.len() % self.decimate, 0);

        // Allocate output
        let mut y: Vec<i16> = Vec::with_capacity(x.len() / self.decimate);
        unsafe { y.set_len(x.len() / self.decimate) };

        // Grab pointers to various things
        let mut delay_idx = self.delay_idx;
        let delay_len = self.delay.len() as isize;
        let delay_p = &mut self.delay[0] as *mut i32;
        let out_p = &mut y[0] as *mut i16;
        let mut in_p = &x[0] as *const i16;
        let ylen = y.len() as isize;
        let decimate = self.decimate as isize;
        let tap0 = &self.taps[0] as *const i32;

        // Process each actually generated output sample
        for k in 0..ylen {

            // Feed the delay line, fast-forwarding through the skipped samples
            for _ in 0..decimate {
                unsafe {
                    *delay_p.offset(delay_idx) = *in_p as i32;
                    in_p = in_p.offset(1);
                    delay_idx -= 1;
                    if delay_idx == -1 {
                        delay_idx = delay_len - 1;
                    }
                }
            }

            // Compute the multiply-accumulate for actual samples.
            let mut acc: i32 = 0;
            let mut tap_p = tap0;

            // First the index to the end of the buffer
            for idx in (delay_idx + 1)..delay_len {
                unsafe {
                    acc += *tap_p * *delay_p.offset(idx);
                    tap_p = tap_p.offset(1);
                }
            }

            // Then the start to the index of the buffer
            for idx in 0..(delay_idx + 1) {
                unsafe {
                    acc += *tap_p * *delay_p.offset(idx);
                    tap_p = tap_p.offset(1);
                }
            }

            // Save the result, accounting for filter gain
            unsafe { *out_p.offset(k) = (acc >> 15) as i16 };
        }

        // Update index for next time
        self.delay_idx = delay_idx;

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

/// Quantise FIR taps to i16 taps that sum to 32768
pub fn quantise_taps_i16(taps: &Vec<f64>) -> Vec<i16> {
    let sum: f64 = taps.iter().fold(0.0_f64, |acc, &x| acc + x);
    taps.iter().map(|t| (t * 32768.0_f64 / sum) as i16).collect()
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
    use super::{FIR, firwin2, quantise_taps_i16, hamming, dft, idft, irdft};

    #[test]
    fn test_fir_impulse() {
        let taps: Vec<i16> = vec!{8192, 16384, 8192};
        let mut fir = FIR::new(&taps, 1);
        let x: Vec<i16> = vec!{4, 0, 0, 0, 0, 0, 0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{1, 2, 1, 0, 0, 0, 0});
    }

    #[test]
    fn test_fir_decimate() {
        let taps: Vec<i16> = vec!{8192, 8192, 8192, 8192};
        let mut fir = FIR::new(&taps, 2);
        let x: Vec<i16> = vec!{4, 4, 4, 4, 8, 8, 8, 8};
        let y = fir.process(&x);
        assert_eq!(y, vec!{2, 4, 6, 8});
    }

    #[test]
    fn test_fir_continuity() {
        let taps: Vec<i16> = vec!{8192, 16384, 8192};
        let mut fir = FIR::new(&taps, 1);
        let x: Vec<i16> = vec!{4, 0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{1, 2});
        let x: Vec<i16> = vec!{0, 0, 0, 0, 0};
        let y = fir.process(&x);
        assert_eq!(y, vec!{1, 0, 0, 0, 0});
    }

    #[test]
    fn test_fir_compensate() {
        let fir = FIR::cic_compensator(63, 5, 8, 2);
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
    fn test_quantise_taps_i16() {
        let taps = vec!{0.25, 0.5, 0.25};
        assert_eq!(quantise_taps_i16(&taps), vec!{8192, 16384, 8192});
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
