extern crate time;
extern crate sdr;

use sdr::downconverter::Downconverter;

#[cfg_attr(test, allow(dead_code))]
fn main() {
    let n_samples = 1<<23 as usize;
    let n_repeats = 10u64;
    let mut dc = Downconverter::new(2);
    let mut x: Vec<u16> = Vec::with_capacity(n_samples);
    for _ in 0..(n_samples/8) {
        x.push(2048+3); x.push(2048+7); x.push(2048+4); x.push(2048+6);
        x.push(2048-3); x.push(2048-7); x.push(2048-4); x.push(2048-6);
    }

    let t0 = time::precise_time_ns();
    for _ in 0..n_repeats {
        dc.process(&x);
    }
    let t1 = time::precise_time_ns();
    let total_samples = n_samples as f64 * n_repeats as f64;
    let total_time = (t1 - t0) as f64 / 1e9;
    let throughput = total_samples / total_time;
    println!("{} blocks of {} samples, {:.2}Msps",
             n_repeats, n_samples, throughput / 1000000.0_f64);
}
