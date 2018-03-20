use num::{Complex, Num, Zero};

pub type IQ<T> = Complex<T>;

/// A trait for types the FM Demodulator can process.
/// Always IQ<T> in and <T> out.
pub trait SampleType: Copy + Num + Zero {
    type Reg: Copy + Num + Zero;
    unsafe fn set_reg(x: *const Complex<Self>, reg: *mut Complex<Self::Reg>);
    unsafe fn set_out(d0re: *const Self::Reg, d0im: *const Self::Reg,
                      d1re: *const Self::Reg, d1im: *const Self::Reg,
                      d2re: *const Self::Reg, d2im: *const Self::Reg,
                      op: *mut Self);
}

macro_rules! impl_fixed_sampletype {
    ($t: ty, $tt: ty, $bits: expr) => {
        impl SampleType for $t {
            type Reg = $tt;

            #[inline]
            unsafe fn set_reg(x: *const Complex<$t>, reg: *mut Complex<$tt>) {
                (*reg).re = (*x).re as $tt;
                (*reg).im = (*x).im as $tt;
            }

            #[inline]
            unsafe fn set_out(d0re: *const Self::Reg, d0im: *const Self::Reg,
                              d1re: *const Self::Reg, d1im: *const Self::Reg,
                              d2re: *const Self::Reg, d2im: *const Self::Reg,
                              op: *mut Self) {
                let top = *d1re * (*d0im - *d2im) - *d1im * (*d0re - *d2re);
                let bottom = (*d1re * *d1re + *d1im * *d1im) >> $bits;
                if bottom == 0 {
                    *op = top as $t;
                } else {
                    *op = (top / bottom) as $t;
                }
            }
        }
    }
}

macro_rules! impl_floating_sampletype {
    ($t: ty) => {
        impl SampleType for $t {
            type Reg = $t;
            #[inline]
            unsafe fn set_reg(x: *const Complex<$t>, reg: *mut Complex<$t>) {
                *reg = *x;
            }

            #[inline]
            unsafe fn set_out(d0re: *const Self::Reg, d0im: *const Self::Reg,
                              d1re: *const Self::Reg, d1im: *const Self::Reg,
                              d2re: *const Self::Reg, d2im: *const Self::Reg,
                              op: *mut Self) {
                let top = *d1re * (*d0im - *d2im) - *d1im * (*d0re - *d2re);
                let bottom = *d1re * *d1re + *d1im * *d1im;
                if bottom == 0.0 {
                    *op = top;
                } else {
                    *op = top / bottom;
                }
            }
        }
    }
}

impl_fixed_sampletype!(i8, i16, 8-2);
impl_fixed_sampletype!(i16, i32, 16-2);
impl_fixed_sampletype!(i32, i64, 32-2);
impl_floating_sampletype!(f32);
impl_floating_sampletype!(f64);

/// FM signal demodulator
#[derive(Default)]
pub struct FMDemod<T: SampleType> {
    d: [IQ<T::Reg>; 3]
}

impl <T:SampleType> FMDemod<T> {
    /// Create a new FM demodulator
    pub fn new() -> FMDemod<T> {
        let d: [IQ<T::Reg>; 3] = [IQ::new(T::Reg::zero(), T::Reg::zero()); 3];
        FMDemod{d }
    }

    /// FM demodulate input block x, containing baseband IQ data that has been
    /// frequency modulated. Outputs a block of real-valued samples.
    pub fn process(&mut self, x: &[IQ<T>]) -> Vec<T> {
        let n = x.len();
        let mut d = self.d;

        // Allocate output
        let mut out: Vec<T> = Vec::with_capacity(n);
        unsafe { out.set_len(n) };

        // Grab pointers to everything for speedy access
        let xp = &x[0] as *const IQ<T>;
        let op = &mut out[0] as *mut T;
        let d0 = &mut d[0] as *mut IQ<T::Reg>;
        let d1 = &mut d[1] as *mut IQ<T::Reg>;
        let d2 = &mut d[2] as *mut IQ<T::Reg>;
        let d0re = &d[0].re as *const T::Reg;
        let d1re = &d[1].re as *const T::Reg;
        let d2re = &d[2].re as *const T::Reg;
        let d0im = &d[0].im as *const T::Reg;
        let d1im = &d[1].im as *const T::Reg;
        let d2im = &d[2].im as *const T::Reg;

        // FM demodulate. Instantaneous phase is approximated by
        // a three-long differentiator doing a multiplication with the
        // old conjugate signal, around the point of interest.
        //
        // Normalised by incoming signal magnitude, but this can go wrong
        // for zero valued signals so we also add a "tiny" delta.
        for i in 0isize..(n as isize) {
            unsafe {
                *d2 = *d1;
                *d1 = *d0;
                T::set_reg(xp.offset(i), d0);
                T::set_out(d0re, d0im, d1re, d1im, d2re, d2im, op.offset(i));
            }
        }

        // Save differentiator state for next time
        self.d = d;

        out
    }
}

#[test]
fn test_fm_demod() {
    let mut x: Vec<IQ<f32>> = Vec::with_capacity(100);
    for t in 0..25 {
        x.push(IQ::new(
            (0.1*(t as f32)).cos(),
            (0.1*(t as f32)).sin()
        ));
    }
    for t in 25..50 {
        x.push(IQ::new(
            (0.3*(t as f32)).cos(),
            (0.3*(t as f32)).sin()
        ));
    }
    for t in 50..75 {
        x.push(IQ::new(
            (0.2*(t as f32)).cos(),
            (0.2*(t as f32)).sin()
        ));
    }
    for t in 75..100 {
        x.push(IQ::new(
            (0.4*(t as f32)).cos(),
            (0.4*(t as f32)).sin()
        ));
    }
    let mut demod = FMDemod::new();
    let y = demod.process(&x);
    let z: Vec<i32> = y.iter().map(|y| (y * 10.0) as i32).collect();
    assert_eq!(z, vec!{
        0, 0,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        -8, -6,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        12, 11,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        6, 8,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7});
}
