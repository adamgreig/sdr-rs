use std::mem::transmute;
use super::IQ;

/// Downconvert by Fs/4
///
/// Takes an input stream of 12-bit u16 reals and outputs i16 IQ values,
/// with Fs/4 downconverted to DC. Does not decimate or antialias, so the
/// IQ sample rate equals the input sample rate, and you will likely want to
/// antialias and then decimate by at least 2 afterwards.
///
/// Input block length must be a multiple of 4 for this process.
#[derive(Default)]
pub struct RealFs4Downconverter;
impl RealFs4Downconverter {
    pub fn new() -> RealFs4Downconverter {
        RealFs4Downconverter
    }

    pub fn process(&self, x: &[u16]) -> Vec<IQ<i16>> {
        let n = x.len();
        assert_eq!(n % 4, 0);

        // Cast input to i16 (we'll clear the sign bit in the loop, in case)
        let x: &[i16] = unsafe { transmute(x) };

        // Allocate output
        let mut y: Vec<IQ<i16>> = Vec::with_capacity(n);
        unsafe { y.set_len(n) };

        // Grab pointers
        let ip = &x[0] as *const i16;
        let op = &mut y[0] as *mut IQ<i16>;

        // Mask and scale to zero-centered i16 values,
        // and flip alternate signs to perform digital downconversion by Fs/4.
        let mut k = 0isize;
        while k < (n as isize) {
            unsafe {
                (*op.offset(k+0)).re =  (*ip.offset(k+0) & 0xFFF) - 2048;
                (*op.offset(k+0)).im = 0;

                (*op.offset(k+1)).re = 0;
                (*op.offset(k+1)).im = -(*ip.offset(k+1) & 0xFFF) + 2048;
                
                (*op.offset(k+2)).re = -(*ip.offset(k+2) & 0xFFF) + 2048;
                (*op.offset(k+2)).im = 0;
                
                (*op.offset(k+3)).re = 0;
                (*op.offset(k+3)).im =  (*ip.offset(k+3) & 0xFFF) - 2048;
            }
            k += 4;
        }

        y
    }
}

#[cfg(test)]
mod tests {
    use super::super::IQ;
    use super::RealFs4Downconverter;

    #[test]
    fn test_downconvert_fs_4() {
        let dc = RealFs4Downconverter::new();
        let x = vec!{2049u16, 2047, 2050, 2046, 2051, 2045, 2052, 2044};
        let y = dc.process(&x);
        assert_eq!(y, vec!{
            IQ::new(1, 0), IQ::new(0, 1), IQ::new(-2, 0), IQ::new(0, -2),
            IQ::new(3, 0), IQ::new(0, 3), IQ::new(-4, 0), IQ::new(0, -4)
        });
    }
}
