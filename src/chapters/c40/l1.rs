//! Chapter 40, Lesson 1: 3D Engagement
//!
//! Three-dimensional missile-target engagement with weaving target.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rtm1: Vec<f64>,
    pub rtm2: Vec<f64>,
    pub rtm3: Vec<f64>,
    pub rtm: Vec<f64>,
}

/// Run the C40L1 simulation
pub fn run() -> Results {
    let _qpn: i32 = 1; // Navigation mode: 1=proportional navigation
    let tau: f64 = 1.0;
    let w: f64 = 3.0;
    let at: f64 = 193.2;
    let vt: f64 = 1000.0;
    let vm: f64 = 3000.0;
    let xnp: f64 = 3.0;
    let xnclim: f64 = 9999999999.0;

    let mut array_t = Vec::new();
    let mut array_rtm1 = Vec::new();
    let mut array_rtm2 = Vec::new();
    let mut array_rtm3 = Vec::new();
    let mut array_rtm = Vec::new();

    let mut rt3ic = 40000.0;
    while rt3ic >= 500.0 {
        let mut rm1: f64 = 0.0;
        let mut rm2: f64 = 10000.0;
        let mut rm3: f64 = 0.0;
        let mut rt1: f64 = 0.0;
        let mut rt2: f64 = 10000.0;
        let mut rt3 = rt3ic;
        let mut vt1 = -at / w;
        let mut vt2: f64 = 0.0;
        let mut vt3 = -vt;
        let mut t: f64 = 0.0;
        let mut _s: f64 = 0.0;
        let mut h: f64;

        let mut vm1: f64 = 0.0;
        let mut vm2: f64 = 0.0;
        let mut vm3 = vm;
        let mut am1: f64 = 0.0;
        let mut am2: f64 = 0.0;
        let mut am3: f64 = 0.0;

        let mut rtm1 = rt1 - rm1;
        let mut rtm2 = rt2 - rm2;
        let mut rtm3 = rt3 - rm3;
        let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        let vtm3 = vt3 - vm3;
        let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;

        while vc >= 0.0 {
            if rtm < 1000.0 {
                h = 0.0002;
            } else {
                h = 0.01;
            }

            let rt1old = rt1;
            let rt2old = rt2;
            let rt3old = rt3;
            let rm1old = rm1;
            let rm2old = rm2;
            let rm3old = rm3;
            let vm1old = vm1;
            let vm2old = vm2;
            let vm3old = vm3;
            let vt1old = vt1;
            let vt2old = vt2;
            let vt3old = vt3;
            let am1old = am1;
            let am2old = am2;
            let am3old = am3;

            // First derivative evaluation
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm3 = rt3 - rm3;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
            let vtm1 = vt1 - vm1;
            let vtm2 = vt2 - vm2;
            let vtm3 = vt3 - vm3;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;
            let tgo = rtm / vc;
            let at1 = at * (w * t).sin();
            let at2 = at * (w * t).cos();
            let at3: f64 = 0.0;

            // ZEM (Zero Effort Miss) calculation - qpn == 1 always for this simulation
            let (zem1, zem2, zem3) = (rtm1 + vtm1 * tgo, rtm2 + vtm2 * tgo, rtm3 + vtm3 * tgo);

            let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2 + zem3 * rtm3) / rtm;
            let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
            let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;
            let zemper3 = zem3 - zemdotrtm * rtm3 / rtm;
            let am1p = xnp * zemper1 / (tgo * tgo);
            let am2p = xnp * zemper2 / (tgo * tgo);
            let am3p = xnp * zemper3 / (tgo * tgo);
            let am1d = (am1p - am1) / tau;
            let am2d = (am2p - am2) / tau;
            let am3d = (am3p - am3) / tau;

            // Euler step
            rt1 += h * vt1;
            rt2 += h * vt2;
            rt3 += h * vt3;
            rm1 += h * vm1;
            rm2 += h * vm2;
            rm3 += h * vm3;
            vm1 += h * am1;
            vm2 += h * am2;
            vm3 += h * am3;
            vt1 += h * at1;
            vt2 += h * at2;
            vt3 += h * at3;
            am1 += h * am1d;
            am2 += h * am2d;
            am3 += h * am3d;
            t += h;

            // RK2 averaging
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm3 = rt3 - rm3;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
            let vtm1 = vt1 - vm1;
            let vtm2 = vt2 - vm2;
            let vtm3 = vt3 - vm3;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;
            let tgo = if vc > 0.0 { rtm / vc } else { 0.001 };
            let at1 = at * (w * t).sin();
            let at2 = at * (w * t).cos();
            let at3: f64 = 0.0;

            // ZEM (Zero Effort Miss) calculation - qpn == 1 always for this simulation
            let (zem1, zem2, zem3) = (rtm1 + vtm1 * tgo, rtm2 + vtm2 * tgo, rtm3 + vtm3 * tgo);

            let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2 + zem3 * rtm3) / rtm;
            let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
            let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;
            let zemper3 = zem3 - zemdotrtm * rtm3 / rtm;
            let am1p = xnp * zemper1 / (tgo * tgo);
            let am2p = xnp * zemper2 / (tgo * tgo);
            let am3p = xnp * zemper3 / (tgo * tgo);
            let am1d = (am1p - am1) / tau;
            let am2d = (am2p - am2) / tau;
            let am3d = (am3p - am3) / tau;

            rt1 = 0.5 * (rt1old + rt1 + h * vt1);
            rt2 = 0.5 * (rt2old + rt2 + h * vt2);
            rt3 = 0.5 * (rt3old + rt3 + h * vt3);
            rm1 = 0.5 * (rm1old + rm1 + h * vm1);
            rm2 = 0.5 * (rm2old + rm2 + h * vm2);
            rm3 = 0.5 * (rm3old + rm3 + h * vm3);
            vm1 = 0.5 * (vm1old + vm1 + h * am1);
            vm2 = 0.5 * (vm2old + vm2 + h * am2);
            vm3 = 0.5 * (vm3old + vm3 + h * am3);
            vt1 = 0.5 * (vt1old + vt1 + h * at1);
            vt2 = 0.5 * (vt2old + vt2 + h * at2);
            vt3 = 0.5 * (vt3old + vt3 + h * at3);
            am1 = 0.5 * (am1old + am1 + h * am1d);
            am2 = 0.5 * (am2old + am2 + h * am2d);
            am3 = 0.5 * (am3old + am3 + h * am3d);

            // Clamp accelerations
            am1 = am1.clamp(-xnclim, xnclim);
            am2 = am2.clamp(-xnclim, xnclim);
            am3 = am3.clamp(-xnclim, xnclim);

            _s += h;
        }

        array_t.push(t);
        array_rtm1.push(rtm1);
        array_rtm2.push(rtm2);
        array_rtm3.push(rtm3);
        array_rtm.push(rtm);

        rt3ic -= 500.0;
    }

    Results {
        time: array_t,
        rtm1: array_rtm1,
        rtm2: array_rtm2,
        rtm3: array_rtm3,
        rtm: array_rtm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c40l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rtm1.clone(),
        results.rtm2.clone(),
        results.rtm3.clone(),
        results.rtm.clone(),
    ])?;

    let plot_file = format!("{}/c40l1_miss.png", output_dir);
    let config = PlotConfig::new("3D Engagement - Miss Distance")
        .with_labels("Time (s)", "RTM (ft)");

    let series = vec![
        Series::new(results.time.clone(), results.rtm.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C40L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c40l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
