use num::{Complex, Num, Zero};

pub type IQ<T> = Complex<T>;

/// A trait for types the FM Demodulator can process.
/// Always IQ<T> in and <T> out.
pub trait SampleType: Copy + Num + Zero {
    fn small() -> Self;
}
impl SampleType for i8  { fn small() -> i8  { 1 }}
impl SampleType for i16 { fn small() -> i16 { 1 }}
impl SampleType for i32 { fn small() -> i32 { 1 }}
impl SampleType for f32 { fn small() -> f32 { 1e-6 }}
impl SampleType for f64 { fn small() -> f64 { 1e-6 }}

/// FM signal demodulator
pub struct FMDemod<T: SampleType> {
    d: [IQ<T>; 3]
}

impl <T:SampleType> FMDemod<T> {
    /// Create a new FM demodulator
    pub fn new() -> FMDemod<T> {
        let d: [IQ<T>; 3] = [IQ::new(T::zero(), T::zero()); 3];
        FMDemod{d: d}
    }

    /// FM demodulate input block x, containing baseband IQ data that has been
    /// frequency modulated. Outputs a block of real-valued samples.
    pub fn process(&mut self, x: &Vec<IQ<T>>) -> Vec<T> {
        let n = x.len();
        let mut d = self.d;

        // Allocate output
        let mut out: Vec<T> = Vec::with_capacity(n);
        unsafe { out.set_len(n) };

        // Grab pointers to everything for speedy access
        let xp = &x[0] as *const IQ<T>;
        let op = &mut out[0] as *mut T;
        let d0 = &mut d[0] as *mut IQ<T>;
        let d1 = &mut d[1] as *mut IQ<T>;
        let d2 = &mut d[2] as *mut IQ<T>;
        let d0re = &d[0].re as *const T;
        let d1re = &d[1].re as *const T;
        let d2re = &d[2].re as *const T;
        let d0im = &d[0].im as *const T;
        let d1im = &d[1].im as *const T;
        let d2im = &d[2].im as *const T;
        let small = T::small();

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
                *d0 = *xp.offset(i);
                *op.offset(i) = (*d1re * (*d0im - *d2im) - *d1im * (*d0re - *d2re))
                                / (*d1re * *d1re + *d1im * *d1im + small);
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
