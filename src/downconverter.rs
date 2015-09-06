use std::mem::transmute;
pub use super::{Real, IQ};
use cic::CIC;
use fir::FIR;

/// Downconverter
pub struct Downconverter {
    n_cics: usize,
    i_cics: Vec<CIC>,
    q_cics: Vec<CIC>,
    i_firs: Vec<FIR>,
    q_firs: Vec<FIR>,
}

impl Downconverter {
    /// Create a new downconverter that has n CIC filters,
    /// each decimates by 16 (/8 from the CIC and /2 from subsequent FIR).
    pub fn new(n_cics: usize) -> Downconverter {
        let mut i_cics: Vec<CIC> = Vec::with_capacity(n_cics);
        let mut q_cics: Vec<CIC> = Vec::with_capacity(n_cics);
        let mut i_firs: Vec<FIR> = Vec::with_capacity(n_cics);
        let mut q_firs: Vec<FIR> = Vec::with_capacity(n_cics);
        for _ in 0..n_cics {
            i_cics.push(CIC::new(5, 8));
            q_cics.push(CIC::new(5, 8));
            i_firs.push(FIR::cic_compensator(64, 5, 8, 2));
            q_firs.push(FIR::cic_compensator(64, 5, 8, 2));
        }
        Downconverter{
            n_cics: n_cics,
            i_cics: i_cics, q_cics: q_cics,
            i_firs: i_firs, q_firs: q_firs,
        }
    }

    /// Downconvert a block of 12 bit unsigned IF=Fs/4 input to a block of
    /// 16 bit signed baseband IQ samples.
    pub fn process(&mut self, x: &Vec<Real<u16>>) -> Vec<IQ<i16>> {
        // Check we were initialised correctly.
        assert!(self.n_cics > 0);
        assert_eq!(self.i_cics.len(), self.n_cics);
        assert_eq!(self.q_cics.len(), self.n_cics);
        assert_eq!(self.i_firs.len(), self.n_cics);
        assert_eq!(self.q_firs.len(), self.n_cics);

        // IQ mix IF down to 0Hz
        let (mut i, mut q) = Downconverter::convert_fs_4(x);

        // Filter and decimate down to baseband
        for idx in 0..self.n_cics {
            i = self.i_cics[idx].process(&i);
            i = self.i_firs[idx].process(&i);
            q = self.q_cics[idx].process(&q);
            q = self.q_firs[idx].process(&q);
        }

        // Return the combined IQ vector
        Downconverter::combine_i_q(&i, &q)
    }

    /// Convert an input stream of 12-bit u16s to an output of i16 I and Q,
    /// downconverted by Fs/4. Does not decimate or antialias, so you will
    /// want to antialias and then decimate by at least 2 after this.
    fn convert_fs_4(x: &Vec<Real<u16>>) -> (Vec<Real<i16>>, Vec<Real<i16>>) {
        // Block length must be a multiple of 4 for downconversion.
        assert_eq!(x.len() % 4, 0);

        // Store new I and Q arrays separately to start with.
        let mut i: Vec<Real<i16>> = Vec::with_capacity(x.len());
        let mut q: Vec<Real<i16>> = Vec::with_capacity(x.len());
        unsafe { i.set_len(x.len()) };
        unsafe { q.set_len(x.len()) };

        // Cast input to i16s (we'll clear the sign bit in the loop anyway).
        let x: &Vec<Real<i16>> = unsafe { transmute(x) };

        // Mask and scale to get zero-centered i16 (input is 12 bit),
        // and flip signs to perform digital downconversion by Fs/4
        // (hence processing in groups of 4, the period of the LOs).
        // TODO: Use range.step_by when it lands.
        let mut k: usize = 0;
        while k < x.len() {
            unsafe {
                *i.get_unchecked_mut(k  ) =   (x.get_unchecked(k  ) & 0xFFF)
                                            - 2048;
                *q.get_unchecked_mut(k  ) = 0;

                *i.get_unchecked_mut(k+1) = 0;
                *q.get_unchecked_mut(k+1) = -((x.get_unchecked(k+1) & 0xFFF)
                                            - 2048);

                *i.get_unchecked_mut(k+2) = -((x.get_unchecked(k+2) & 0xFFF)
                                            - 2048);
                *q.get_unchecked_mut(k+2) = 0;

                *i.get_unchecked_mut(k+3) = 0;
                *q.get_unchecked_mut(k+3) =   (x.get_unchecked(k+3) & 0xFFF)
                                            - 2048;
            }
            k += 4;
        }

        (i, q)
    }

    /// Combine separate I and Q vectors into a single Vec<IQ>.
    fn combine_i_q(i: &Vec<Real<i16>>, q: &Vec<Real<i16>>) -> Vec<IQ<i16>> {
        // Check I and Q inputs are the same length
        assert_eq!(i.len(), q.len());

        // Combine
        i.iter().zip(q).map(|(i, q)| IQ::new(*i, *q)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::Downconverter;
    use super::IQ;

    #[test]
    fn test_convert_fs_4() {
        let x = vec!{2049u16, 2047, 2050, 2046, 2051, 2045, 2052, 2044};
        let (i, q) = Downconverter::convert_fs_4(&x);
        assert_eq!(i, vec!{1, 0, -2, 0, 3, 0, -4, 0});
        assert_eq!(q, vec!{0, 1, 0, -2, 0, 3, 0, -4});
    }

    #[test]
    fn test_combining_i_q() {
        let i = vec!{5i16, -5};
        let q = vec!{-3i16, 3};
        let iq = Downconverter::combine_i_q(&i, &q);
        assert_eq!(iq, vec!{IQ::new(5i16, -3i16), IQ::new(-5i16, 3i16)});
    }
}
