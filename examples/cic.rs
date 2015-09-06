extern crate time;
extern crate sdr;
use sdr::CIC;
use sdr::IQ;

#[cfg_attr(test, allow(dead_code))]
fn main() {
    bench_scalar();
    bench_complex();
}

#[cfg_attr(test, allow(dead_code))]
fn bench_scalar() {
    let q = 5usize;
    let r = 16usize;
    let n_samples = 1<<24 as usize;
    let n_repeats = 10u64;

    let mut cic = CIC::new(q, r, 12);
    let mut x: Vec<i16> = Vec::with_capacity(n_samples);
    for _ in 0..(n_samples/8) {
        x.push(3); x.push(7); x.push(4); x.push(6);
        x.push(3); x.push(7); x.push(4); x.push(6);
    }
    let t0 = time::precise_time_ns();
    for _ in 0..n_repeats {
        cic.process(&x);
    }
    let t1 = time::precise_time_ns();
    let total_samples = n_samples as f64 * n_repeats as f64;
    let total_time = (t1 - t0) as f64 / 1e9;
    let throughput = total_samples / total_time;
    println!("i16:     {} blocks of {} samples, {:.2}Msps",
             n_repeats, n_samples, throughput / 1000000.0_f64);
}

#[cfg_attr(test, allow(dead_code))]
fn bench_complex() {
    let q = 5usize;
    let r = 16usize;
    let n_samples = 1<<24 as usize;
    let n_repeats = 10u64;

    let mut cic = CIC::new(q, r, 12);
    let mut x: Vec<IQ<i16>> = Vec::with_capacity(n_samples);
    for _ in 0..(n_samples/8) {
        x.push(IQ{re:3, im:-3}); x.push(IQ{re:7, im:-7});
        x.push(IQ{re:4, im:-4}); x.push(IQ{re:6, im:-6});
        x.push(IQ{re:3, im:-3}); x.push(IQ{re:7, im:-7});
        x.push(IQ{re:4, im:-4}); x.push(IQ{re:6, im:-6});
    }
    let t0 = time::precise_time_ns();
    for _ in 0..n_repeats {
        cic.process(&x);
    }
    let t1 = time::precise_time_ns();
    let total_samples = n_samples as f64 * n_repeats as f64;
    let total_time = (t1 - t0) as f64 / 1e9;
    let throughput = total_samples / total_time;
    println!("IQ<i16>: {} blocks of {} samples, {:.2}Msps",
             n_repeats, n_samples, throughput / 1000000.0_f64);
}
