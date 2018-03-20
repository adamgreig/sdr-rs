use num::{Zero, Complex};

/// Implement the underlying operations and types required by allowable CIC
/// sample types. Must be some form of wrapping arithmetic.
pub trait SampleType: Copy {
    type RegType: Zero + Copy;
    fn check_width(q: usize, r: usize, bw: usize);
    fn add(a: Self::RegType, b: Self::RegType) -> Self::RegType;
    fn sub(a: Self::RegType, b: Self::RegType) -> Self::RegType;
    fn to_reg(x: Self) -> Self::RegType;
    fn from_out(out: Self::RegType, gain_shift: usize) -> Self;
}

/// Implement `SampleType` for scalars like i16.
/// $t is the sample type, $tt the register type,
/// $bw_t the bit width of type t, and $bw_tt the bit width of $tt.
macro_rules! impl_scalar_sampletype {
    ($t:ty, $tt:ty) => {
        impl SampleType for $t {
            type RegType = $tt;
            #[inline] fn check_width(q: usize, r: usize, bw: usize) {
                let bw = bw as f32;
                let bw_t = <$t>::zero().count_zeros() as f32;
                let bw_tt = <$tt>::zero().count_zeros() as f32;
                assert!(bw <= bw_t);
                assert!(bw + ((q as f32) * (r as f32).log2()).ceil() <= bw_tt);
            }
            #[inline] fn add(a: $tt, b: $tt) -> $tt {
                a.wrapping_add(b)
            }
            #[inline] fn sub(a: $tt, b: $tt) -> $tt {
                a.wrapping_sub(b)
            }
            #[inline] fn to_reg(x: $t) -> $tt {
                x as $tt
            }
            #[inline] fn from_out(out: $tt, gain_shift: usize) -> $t {
                (out >> gain_shift) as $t
            }
        }
    }
}

/// Implement `SampleType` for complex samples like `Complex<i16>`.
/// $t is the sample type, $tt the register type,
/// $bw_t the bit width of type t, and $bw_tt the bit width of $tt.
macro_rules! impl_complex_sampletype {
    ($t:ty, $tt:ty) => {
        impl SampleType for Complex<$t> {
            type RegType = Complex<$tt>;
            #[inline] fn check_width(q: usize, r: usize, bw: usize) {
                let bw = bw as f32;
                let bw_t = <$t>::zero().count_zeros() as f32;
                let bw_tt = <$tt>::zero().count_zeros() as f32;
                assert!(bw <= bw_t);
                assert!(bw + ((q as f32) * (r as f32).log2()).ceil() <= bw_tt);
            }
            #[inline]
            fn add(a: Complex<$tt>, b: Complex<$tt>) -> Complex<$tt> {
                Complex{ re: a.re.wrapping_add(b.re),
                         im: a.im.wrapping_add(b.im) }
            }
            #[inline]
            fn sub(a: Complex<$tt>, b: Complex<$tt>) -> Complex<$tt> {
                Complex{ re: a.re.wrapping_sub(b.re),
                         im: a.im.wrapping_sub(b.im) }
            }
            #[inline] fn to_reg(x: Complex<$t>) -> Complex<$tt> {
                Complex{ re: x.re as $tt, im: x.im as $tt }
            }
            #[inline] fn from_out(out: Complex<$tt>, shift: usize)
                -> Complex<$t>
            {
                Complex{ re: (out.re >> shift) as $t,
                         im: (out.im >> shift) as $t }
            }
        }
    }
}

macro_rules! impl_sampletype {
    ($t:ty, $tt:ty) => {
        impl_scalar_sampletype!($t, $tt);
        impl_complex_sampletype!($t, $tt);
    }
}

impl_sampletype!(i8, i16);
impl_sampletype!(i16, i32);
impl_sampletype!(i32, i64);

/// Qth order CIC filter with decimate-by-R (and D=R).
pub struct CIC<T: SampleType> {
    q: usize,
    r: usize,
    intg: Vec<T::RegType>,
    comb: Vec<T::RegType>,
    gain_shift: usize
}

impl <T: SampleType> CIC<T> {
    /// Create a new CIC filter of order q and downsampling ratio r.
    ///
    /// Bitwidth is the number of bits actually used by your input data,
    /// so must be equal to or less than the size of the data type.
    ///
    /// Note that bitwidth + ceil(Q log2 R) must <= 32.
    pub fn new(q: usize, r: usize, bw: usize) -> CIC<T> {
        // Verify that bw + ceil(Q log2 R) <= reg width.
        T::check_width(q, r, bw);

        // Allocate the integrator and comb registers.
        let mut intg = Vec::with_capacity(q);
        let mut comb = Vec::with_capacity(q+1);

        for _ in 0..q {
            intg.push(T::RegType::zero());
            comb.push(T::RegType::zero());
        }

        // One final comb register for holding the output.
        comb.push(T::RegType::zero());

        // Compute the filter gain, Q**R, as an equivalent bit shift
        let gain_shift: usize = (q as f32 * (r as f32).log2()).ceil() as usize;

        CIC { q, r, intg, comb, gain_shift }
    }

    /// Run the CIC filter over a block x,
    /// returning the filtered and decimated output.
    pub fn process(&mut self, x: &[T]) -> Vec<T> {
        // Check we were initialised correctly
        assert!(self.q > 0 && self.r > 0);
        assert_eq!(self.intg.len(), self.q);
        assert_eq!(self.comb.len(), self.q + 1);

        // To decimate by R we need a multiple of R in the input.
        assert_eq!(x.len() % self.r, 0);

        // Output will have 1/R the number of samples of the input.
        let mut y: Vec<T> = Vec::with_capacity(x.len() / self.r);
        unsafe { y.set_len(x.len() / self.r) };

        // Hold on to some pointers for more dangerous access later.
        let intg_p = &mut self.intg[0] as *mut T::RegType;
        let comb_p = &mut self.comb[0] as *mut T::RegType;
        let out_p = &mut y[0] as *mut T;
        let in_p = &x[0] as *const T;
        let r = self.r as isize;
        let q = self.q as isize;
        let ylen = y.len() as isize;
        let gain_shift = self.gain_shift;

        // Push samples through the chain.
        // It's a fair bit faster to loop over the outputs and have a tighter
        // loop over the inputs.
        for k in 0..ylen {
            for o in 0..r {
                // Add up the integrators. Note that this cascaded approach
                // adds additional time delay (but not much!) but is easier to
                // compute.
                unsafe { 
                    *intg_p =
                        T::add(*intg_p, T::to_reg(*in_p.offset(k*r + o)));

                    for l in 1isize..q {
                        *intg_p.offset(l) =
                            T::add(*intg_p.offset(l), *intg_p.offset(l - 1));
                    }

                }
            }

            // Run the comb section at 1/R the sample rate
            // Each comb register is set to the output of the decimator,
            // minus all the combs before itself
            for l in 0..(q + 1) {
                let l = q - l;
                unsafe {
                    *comb_p.offset(l) = *intg_p.offset(q - 1);

                    for m in 0isize..l {
                        *comb_p.offset(l) =
                            T::sub(*comb_p.offset(l), *comb_p.offset(m));
                    }
                }
            }

            // The output is the final "output" comb register, scaled to
            // give unity filter gain.
            unsafe {
                *out_p.offset(k) = T::from_out(*comb_p.offset(q), gain_shift);
            }
        }

        y
    }
}

#[cfg(test)]
mod tests {
    use super::CIC;
    use num::Complex;

    /// Try 1st order R=2 CIC, a simple moving average filter in practice
    #[test]
    fn test_cic_1_2() {
        let mut cic = CIC::new(1, 2, 12);
        let x = vec!{1i16, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        let y = cic.process(&x);
        assert_eq!(y, vec!{1i16, 1, 1, 1, 1, 1, 1, 1});

        let x = vec!{1i16, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3};
        let mut cic = CIC::new(1, 2, 12);
        let y = cic.process(&x);
        assert_eq!(y, vec!{2i16, 2, 2, 2, 2, 2, 2, 2});
    }

    /// Try a bigger filter
    #[test]
    fn test_cic_2_4() {
        let x = vec!{1i16, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3};
        let mut cic = CIC::new(2, 4, 12);
        let y = cic.process(&x);
        assert_eq!(y, vec!{1i16, 2, 2, 2});
    }

    /// Now we're cooking with gas!
    #[test]
    fn test_cic_5_8() {
        let mut x: Vec<i16> = Vec::with_capacity(8*8);
        for _ in 0..16 {
            x.push(2);
            x.push(4);
            x.push(6);
            x.push(4);
        }
        let mut cic = CIC::new(5, 8, 12);
        let y = cic.process(&x);
        assert_eq!(y, vec!{0i16, 1, 3, 3, 4, 4, 4, 4});
    }

    /// Try a different sample size
    #[test]
    fn test_cic_i32() {
        let x = vec!{1i32, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3};
        let mut cic = CIC::new(2, 4, 5);
        let y = cic.process(&x);
        assert_eq!(y, vec!{1i32, 2, 2, 2});
    }

    /// Try a complex sample
    #[test]
    fn test_cic_complex() {
        let x = vec!{
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
            Complex{ re: 1i16, im: 4i16}, Complex{ re: 3i16, im: 8i16},
        };
        let mut cic = CIC::new(2, 4, 5);
        let y = cic.process(&x);
        assert_eq!(y, vec!{
            Complex{ re: 1i16, im: 3i16 },
            Complex{ re: 2i16, im: 6i16 },
            Complex{ re: 2i16, im: 6i16 },
            Complex{ re: 2i16, im: 6i16 },
        });
    }
}
