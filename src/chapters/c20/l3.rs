//! Chapter 20, Lesson 3: Target Displacement Nonlinear
//!
//! Nonlinear simulation for target displacement at various homing times.

use crate::save_data;

pub struct Results {
    pub thom: Vec<f64>,
    pub rtmp: Vec<f64>,
    pub theory: Vec<f64>,
}

/// Run the C20L3 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let displace: f64 = 200.0;
    let tau: f64 = 1.0;
    let vm: f64 = 3000.0;
    let vt: f64 = 1000.0;

    let mut array_thom = Vec::new();
    let mut array_rtmp = Vec::new();
    let mut array_theory = Vec::new();

    let mut thom = 0.1;
    while thom <= 5.0 + 0.00001 {
        let mut rm1: f64 = 0.0;
        let mut rm2: f64 = 1000.0;
        let mut rt1: f64 = 20000.0;
        let mut rt2: f64 = 1000.0;
        let mut qswitch = false;
        let vt1 = -vt;
        let vt2: f64 = 0.0;

        let mut rtm1 = rt1 - rm1;
        let mut rtm2 = rt2 - rm2;
        let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let mut vm1 = vm;
        let mut vm2: f64 = 0.0;
        let mut vtm1 = vt1 - vm1;
        let mut vtm2 = vt2 - vm2;
        let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let mut xlamh: f64 = 0.0;
        let mut h = 0.001;

        while vc > 0.0 {
            let tgo = rtm / vc;
            if tgo < 0.3 {
                h = 0.00001;
            }

            if tgo <= thom && !qswitch {
                qswitch = true;
                rt2 += displace;
            }

            let rt1old = rt1;
            let rt2old = rt2;
            let rm1old = rm1;
            let rm2old = rm2;
            let vm1old = vm1;
            let vm2old = vm2;
            let xlamhold = xlamh;

            // First Euler step
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            vtm1 = vt1 - vm1;
            vtm2 = vt2 - vm2;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
            let xlam = rtm2.atan2(rtm1);
            let xlamhd = (xlam - xlamh) / tau;
            let xnc = xnp * vc * xlamhd;
            let am1 = -xnc * xlam.sin();
            let am2 = xnc * xlam.cos();

            rt1 += h * vt1;
            rt2 += h * vt2;
            rm1 += h * vm1;
            rm2 += h * vm2;
            vm1 += h * am1;
            vm2 += h * am2;
            xlamh += h * xlamhd;
            // t += h; // time tracking (not used)

            // Second evaluation for RK2
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            vtm1 = vt1 - vm1;
            vtm2 = vt2 - vm2;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
            let xlam = rtm2.atan2(rtm1);
            let xlamhd = (xlam - xlamh) / tau;
            let xnc = xnp * vc * xlamhd;
            let am1 = -xnc * xlam.sin();
            let am2 = xnc * xlam.cos();

            // RK2 averaging
            rt1 = 0.5 * (rt1old + rt1 + h * vt1);
            rt2 = 0.5 * (rt2old + rt2 + h * vt2);
            rm1 = 0.5 * (rm1old + rm1 + h * vm1);
            rm2 = 0.5 * (rm2old + rm2 + h * vm2);
            vm1 = 0.5 * (vm1old + vm1 + h * am1);
            vm2 = 0.5 * (vm2old + vm2 + h * am2);
            xlamh = 0.5 * (xlamhold + xlamh + h * xlamhd);
        }

        let rtmp = if rtm2 > 0.0 { rtm } else { -rtm };
        let x = thom / tau;
        let theory = displace * (-x).exp() * (1.0 - 2.0 * x + 0.5 * x * x);

        array_thom.push(thom);
        array_rtmp.push(rtmp);
        array_theory.push(theory);

        thom += 0.1;
    }

    Results {
        thom: array_thom,
        rtmp: array_rtmp,
        theory: array_theory,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.thom.clone(),
        results.rtmp.clone(),
        results.theory.clone(),
    ];
    let data_file = format!("{}/c20l3_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C20L3: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l3_runs() {
        let results = run();
        assert!(!results.thom.is_empty());
    }
}
