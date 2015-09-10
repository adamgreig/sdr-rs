use num::{Complex, Num, Zero};

pub type IQ<T> = Complex<T>;

pub trait SampleType: Copy + Num + Zero {}
impl SampleType for i8 {}
impl SampleType for i16 {}
impl SampleType for i32 {}
impl SampleType for f32 {}
impl SampleType for f64 {}

pub fn fm_demod<T: SampleType>(x: Vec<IQ<T>>) -> Vec<T> {
    let n = x.len();
    let mut d: [IQ<T>; 3] = [IQ::new(T::zero(), T::zero()); 3];
    let mut out: Vec<T> = Vec::with_capacity(n);
    unsafe { out.set_len(n) };
    
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

    unsafe { *d0 = *xp; *d1 = *xp; }

    for i in 0isize..(n as isize) {
        unsafe {
            *d2 = *d1;
            *d1 = *d0;
            *d0 = *xp.offset(i);
            *op.offset(i) = (*d1re * (*d0im - *d2im) - *d1im * (*d0re - *d2re))
                            / (*d1re * *d1re + *d1im * *d1im);
        }
    }

    out
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
    let y = fm_demod(x);
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
