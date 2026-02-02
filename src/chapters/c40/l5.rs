//! Chapter 40, Lesson 5: 3D Reentry Guidance Simulation
//!
//! 3D engagement with optimal guidance and flight path constraints.

use crate::save_data;
use crate::utils::launch_logic;

pub struct Results {
    pub time: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub rm3k: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub rt3k: Vec<f64>,
    pub gamp_deg: Vec<f64>,
    pub gamfp_deg: Vec<f64>,
    pub gamy_deg: Vec<f64>,
    pub gamfy_deg: Vec<f64>,
}

/// Run the C40L5 simulation
pub fn run() -> Results {
    let xntg: f64 = 4.0;
    let vt: f64 = 1000.0;
    let vm: f64 = 3000.0;
    let mut rm1: f64 = 0.0;
    let mut rm2: f64 = 10000.0;
    let mut rm3: f64 = -1000.0;
    let rt1: f64 = 30000.0;
    let rt2: f64 = 10000.0;
    let rt3: f64 = 0.0;
    let gamfpdeg: f64 = -30.0;
    let gamfydeg: f64 = 20.0;
    let xnt = 32.2 * xntg;
    let mut beta: f64 = 0.0;
    let mut vt1 = -vt * beta.cos();
    let mut vt2 = vt * beta.sin();
    let mut vt3: f64 = 0.0;
    let gamfp = gamfpdeg / 57.3;
    let gamfy = gamfydeg / 57.3;

    // Launch logic to find initial missile velocity
    let (mut vm1, mut vm2, mut vm3, _tf) = launch_logic(rm1, rm2, rm3, rt1, rt2, rt3, vt1, vt2, vt3, vm);

    let mut h: f64;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut rt1 = rt1;
    let mut rt2 = rt2;
    let mut rt3 = rt3;

    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm3 = rt3 - rm3;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
    let mut vtm1 = vt1 - vm1;
    let mut vtm2 = vt2 - vm2;
    let mut vtm3 = vt3 - vm3;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;

    let vm1f = vm * gamfp.cos() * gamfy.cos();
    let vm2f = vm * gamfp.sin();
    let vm3f = vm * gamfp.cos() * gamfy.sin();

    let mut array_t = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_rm3k = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rt3k = Vec::new();
    let mut array_gampdeg = Vec::new();
    let mut array_gamfpdeg = Vec::new();
    let mut array_gamydeg = Vec::new();
    let mut array_gamfydeg = Vec::new();
    let mut array_xncg = Vec::new();

    let mut _n = 0;
    let mut gampdeg: f64 = 0.0;
    let mut gamydeg: f64 = 0.0;

    while vc >= 0.0 {
        if rtm < 1000.0 {
            h = 0.00001;
        } else {
            h = 0.0001;
        }

        let betaold = beta;
        let rt1old = rt1;
        let rt2old = rt2;
        let rt3old = rt3;
        let rm1old = rm1;
        let rm2old = rm2;
        let rm3old = rm3;
        let vm1old = vm1;
        let vm2old = vm2;
        let vm3old = vm3;

        // First derivative evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm3 = rt3 - rm3;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
        vtm1 = vt1 - vm1;
        vtm2 = vt2 - vm2;
        vtm3 = vt3 - vm3;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;
        let tgo = rtm / vc;
        let betad = xnt / vt;
        let xnt1 = xnt * beta.sin();
        let xnt2 = xnt * beta.cos();
        let xnt3: f64 = 0.0;
        vt1 = -vt * beta.cos();
        vt2 = vt * beta.sin();
        vt3 = 0.0;

        let vm_mag = (vm1 * vm1 + vm2 * vm2 + vm3 * vm3).sqrt();
        let xnc1 = (6.0 * rtm1 + 6.0 * vtm1 * tgo) / (tgo * tgo) + 2.0 * (vm1 - vm1f) / tgo + xnt1;
        let xnc2 = (6.0 * rtm2 + 6.0 * vtm2 * tgo) / (tgo * tgo) + 2.0 * (vm2 - vm2f) / tgo + xnt2;
        let xnc3 = (6.0 * rtm3 + 6.0 * vtm3 * tgo) / (tgo * tgo) + 2.0 * (vm3 - vm3f) / tgo + xnt3;
        let xncg = (xnc1 * xnc1 + xnc2 * xnc2 + xnc3 * xnc3).sqrt() / 32.2;
        let xncdotvm = (xnc1 * vm1 + xnc2 * vm2 + xnc3 * vm3) / vm_mag;
        let am1 = xnc1 - xncdotvm * vm1 / vm_mag;
        let am2 = xnc2 - xncdotvm * vm2 / vm_mag;
        let am3 = xnc3 - xncdotvm * vm3 / vm_mag;

        // Euler step
        beta += h * betad;
        rt1 += h * vt1;
        rt2 += h * vt2;
        rt3 += h * vt3;
        rm1 += h * vm1;
        rm2 += h * vm2;
        rm3 += h * vm3;
        vm1 += h * am1;
        vm2 += h * am2;
        vm3 += h * am3;
        t += h;

        // Second evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm3 = rt3 - rm3;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
        vtm1 = vt1 - vm1;
        vtm2 = vt2 - vm2;
        vtm3 = vt3 - vm3;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;
        let tgo = rtm / vc;
        let betad = xnt / vt;
        let xnt1 = xnt * beta.sin();
        let xnt2 = xnt * beta.cos();
        let xnt3: f64 = 0.0;
        vt1 = -vt * beta.cos();
        vt2 = vt * beta.sin();
        vt3 = 0.0;

        let vm_mag = (vm1 * vm1 + vm2 * vm2 + vm3 * vm3).sqrt();
        let xnc1 = (6.0 * rtm1 + 6.0 * vtm1 * tgo) / (tgo * tgo) + 2.0 * (vm1 - vm1f) / tgo + xnt1;
        let xnc2 = (6.0 * rtm2 + 6.0 * vtm2 * tgo) / (tgo * tgo) + 2.0 * (vm2 - vm2f) / tgo + xnt2;
        let xnc3 = (6.0 * rtm3 + 6.0 * vtm3 * tgo) / (tgo * tgo) + 2.0 * (vm3 - vm3f) / tgo + xnt3;
        let xncdotvm = (xnc1 * vm1 + xnc2 * vm2 + xnc3 * vm3) / vm_mag;
        let am1 = xnc1 - xncdotvm * vm1 / vm_mag;
        let am2 = xnc2 - xncdotvm * vm2 / vm_mag;
        let am3 = xnc3 - xncdotvm * vm3 / vm_mag;
        gamydeg = 57.3 * vm3.atan2(vm1);
        gampdeg = 57.3 * vm2.atan2((vm1 * vm1 + vm3 * vm3).sqrt());

        // RK2 averaging
        beta = 0.5 * (betaold + beta + h * betad);
        rt1 = 0.5 * (rt1old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2old + rt2 + h * vt2);
        rt3 = 0.5 * (rt3old + rt3 + h * vt3);
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        rm3 = 0.5 * (rm3old + rm3 + h * vm3);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);
        vm3 = 0.5 * (vm3old + vm3 + h * am3);

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            _n += 1;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;
            let rt3k = rt3 / 1000.0;
            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;
            let rm3k = rm3 / 1000.0;

            array_t.push(t);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_rm3k.push(rm3k);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rt3k.push(rt3k);
            array_gampdeg.push(gampdeg);
            array_gamfpdeg.push(gamfpdeg);
            array_gamydeg.push(gamydeg);
            array_gamfydeg.push(gamfydeg);
            array_xncg.push(xncg);
        }
    }

    // Final point
    let rt1k = rt1 / 1000.0;
    let rt2k = rt2 / 1000.0;
    let rt3k = rt3 / 1000.0;
    let rm1k = rm1 / 1000.0;
    let rm2k = rm2 / 1000.0;
    let rm3k = rm3 / 1000.0;

    array_t.push(t);
    array_rm1k.push(rm1k);
    array_rm2k.push(rm2k);
    array_rm3k.push(rm3k);
    array_rt1k.push(rt1k);
    array_rt2k.push(rt2k);
    array_rt3k.push(rt3k);
    array_gampdeg.push(gampdeg);
    array_gamfpdeg.push(gamfpdeg);
    array_gamydeg.push(gamydeg);
    array_gamfydeg.push(gamfydeg);

    Results {
        time: array_t,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        rm3k: array_rm3k,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        rt3k: array_rt3k,
        gamp_deg: array_gampdeg,
        gamfp_deg: array_gamfpdeg,
        gamy_deg: array_gamydeg,
        gamfy_deg: array_gamfydeg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c40l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.rm3k.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rt3k.clone(),
        results.gamp_deg.clone(),
        results.gamfp_deg.clone(),
        results.gamy_deg.clone(),
        results.gamfy_deg.clone(),
    ])?;

    println!("C40L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c40l5_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
