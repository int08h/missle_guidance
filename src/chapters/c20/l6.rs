//! Chapter 20, Listing 6: 2D Engagement with LOS Filter (Sweep THOM)
//!
//! Simulates 2D engagement with target displacement using a first-order
//! lag filter on the line-of-sight angle. Sweeps through different
//! homing times (THOM). Uses filtered LOS for guidance.

use crate::save_data;

pub struct Results {
    pub thom: Vec<f64>,
    pub rtmp: Vec<f64>,
}

/// Run the C20L6 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let displace: f64 = 200.0;
    let tau: f64 = 1.0;
    let vm: f64 = 3000.0;
    let vt: f64 = 1000.0;

    let mut array_thom = Vec::new();
    let mut array_rtmp = Vec::new();

    // Sweep THOM from 0.1 to 5.0 in steps of 0.1
    let mut thom_val: f64 = 0.1;
    while thom_val <= 5.0 + 0.00001 {
        // Initial conditions
        let mut rm1: f64 = 0.0;
        let mut rm2: f64 = 1000.0;
        let mut rt1: f64 = 20000.0;
        let mut rt2: f64 = 1000.0;
        let mut qswitch = false;
        let vt1 = -vt;
        let vt2: f64 = 0.0;
        let mut _t: f64 = 0.0;

        let mut rtm1 = rt1 - rm1;
        let mut rtm2 = rt2 - rm2;
        let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let mut xlam: f64;

        let mut vm1 = vm;
        let mut vm2: f64 = 0.0;

        let mut vtm1 = vt1 - vm1;
        let mut vtm2 = vt2 - vm2;
        let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let mut tgo: f64;
        let mut d: f64 = 0.0;
        let mut h: f64 = 0.001;

        while vc > 0.0 {
            tgo = rtm / vc;
            if tgo < 0.3 {
                h = 0.00001;
            }

            // Target jink at specified time-to-go
            if tgo <= thom_val && !qswitch {
                qswitch = true;
                rt2 += displace;
                rtm2 = rt2 - rm2;
                xlam = rtm2.atan2(rtm1);
                d = xlam;
            }

            let rt1old = rt1;
            let rt2old = rt2;
            let rm1old = rm1;
            let rm2old = rm2;
            let vm1old = vm1;
            let vm2old = vm2;
            let dold = d;

            // First derivative evaluation
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

            vtm1 = vt1 - vm1;
            vtm2 = vt2 - vm2;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

            xlam = rtm2.atan2(rtm1);
            let _eps = xlam - d;
            let mut dd = (xlam - d) / tau;
            let mut xnc = xnp * vc * dd;
            let mut am1 = -xnc * xlam.sin();
            let mut am2 = xnc * xlam.cos();

            // Euler step
            rt1 += h * vt1;
            rt2 += h * vt2;
            rm1 += h * vm1;
            rm2 += h * vm2;
            vm1 += h * am1;
            vm2 += h * am2;
            d += h * dd;
            _t += h;

            // Second derivative evaluation
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

            vtm1 = vt1 - vm1;
            vtm2 = vt2 - vm2;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

            xlam = rtm2.atan2(rtm1);
            dd = (xlam - d) / tau;
            xnc = xnp * vc * dd;
            am1 = -xnc * xlam.sin();
            am2 = xnc * xlam.cos();

            // RK2 averaging
            rt1 = 0.5 * (rt1old + rt1 + h * vt1);
            rt2 = 0.5 * (rt2old + rt2 + h * vt2);
            rm1 = 0.5 * (rm1old + rm1 + h * vm1);
            rm2 = 0.5 * (rm2old + rm2 + h * vm2);
            vm1 = 0.5 * (vm1old + vm1 + h * am1);
            vm2 = 0.5 * (vm2old + vm2 + h * am2);
            d = 0.5 * (dold + d + h * dd);
        }

        // Compute signed miss distance
        let rtmp = if rtm2 > 0.0 { rtm } else { -rtm };

        array_thom.push(thom_val);
        array_rtmp.push(rtmp);

        thom_val += 0.1;
    }

    Results {
        thom: array_thom,
        rtmp: array_rtmp,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.thom.clone(),
        results.rtmp.clone(),
    ];
    let data_file = format!("{}/c20l6_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C20L6: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l6_runs() {
        let results = run();
        assert!(!results.thom.is_empty());
    }
}
