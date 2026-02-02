//! Chapter 32, Lesson 4: Rolling Airframe Guidance (Closed Loop)
//!
//! Simulates a rolling airframe missile with closed-loop predictive
//! guidance using roll rate control.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xm: Vec<f64>,
    pub ym: Vec<f64>,
    pub xt: Vec<f64>,
    pub yt: Vec<f64>,
    pub phiddeg: Vec<f64>,
}

/// PREDICT1 function - predict final position
#[allow(clippy::too_many_arguments)]
fn predict1(tp: f64, xp: f64, yp: f64, xdp: f64, ydp: f64, thetp: f64, tf: f64, thadp: f64, thbdp: f64, acc: f64, hp: f64, thdmax: f64, accerr: f64) -> (f64, f64) {
    let h = hp * 10.0;
    let mut t = tp;
    let mut x = xp;
    let mut y = yp;
    let mut xd = xdp;
    let mut yd = ydp;
    let mut thet = thetp;
    let thad = thadp;
    let thbd = thbdp;

    while t <= tf - 0.00001 {
        let xold = x;
        let yold = y;
        let xdold = xd;
        let ydold = yd;
        let thetold = thet;

        // First derivative evaluation
        let slope = (thbd - thad) / tf;
        let bint = thbd - slope * tf;
        let mut thetd = slope * t + bint;
        if thetd > thdmax {
            thetd = thdmax;
        }
        if thetd < -thdmax {
            thetd = -thdmax;
        }
        let xdd = (acc + accerr) * thet.cos();
        let ydd = (acc + accerr) * thet.sin();

        // Euler step
        x += h * xd;
        y += h * yd;
        xd += h * xdd;
        yd += h * ydd;
        thet += h * thetd;
        t += h;

        // Second derivative evaluation
        let slope = (thbd - thad) / tf;
        let bint = thbd - slope * tf;
        let mut thetd = slope * t + bint;
        if thetd > thdmax {
            thetd = thdmax;
        }
        if thetd < -thdmax {
            thetd = -thdmax;
        }
        let xdd = (acc + accerr) * thet.cos();
        let ydd = (acc + accerr) * thet.sin();

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * h * xd;
        y = (yold + y) / 2.0 + 0.5 * h * yd;
        xd = (xdold + xd) / 2.0 + 0.5 * h * xdd;
        yd = (ydold + yd) / 2.0 + 0.5 * h * ydd;
        thet = (thetold + thet) / 2.0 + 0.5 * h * thetd;
    }

    (x, y)
}

/// Run the C32L4 simulation
pub fn run() -> Results {
    let gain: f64 = 1.0;
    let xt: f64 = -4000.0;
    let yt: f64 = 5000.0;
    let phiaddeg: f64 = 10.0;
    let phibddeg: f64 = 5.0;
    let phidmaxdeg: f64 = 15.0;
    let tf: f64 = 100.0;
    let accerr: f64 = 0.0;
    let h: f64 = 0.01;
    let ts: f64 = 0.1;
    let xnc: f64 = 12.0;
    let phidmax = phidmaxdeg / 57.3;

    let mut phi: f64 = 0.0;
    let mut xm: f64 = 0.0;
    let mut ym: f64 = 0.0;
    let mut xmd: f64 = 0.0;
    let mut ymd: f64 = 0.0;
    let mut phiad = phiaddeg / 57.3;
    let mut phibd = phibddeg / 57.3;
    let mut s: f64 = 0.0;
    let mut t: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xm = Vec::new();
    let mut array_ym = Vec::new();
    let mut array_xt = Vec::new();
    let mut array_yt = Vec::new();
    let mut array_phiddeg = Vec::new();

    while t <= tf - 0.00001 {
        let xmold = xm;
        let ymold = ym;
        let xmdold = xmd;
        let ymdold = ymd;
        let phiold = phi;

        // First derivative evaluation
        let slope = (phibd - phiad) / tf;
        let bint = phibd - slope * tf;
        let mut phid = slope * t + bint;
        if phid > phidmax {
            phid = phidmax;
        }
        if phid < -phidmax {
            phid = -phidmax;
        }
        let xmdd = xnc * phi.cos();
        let ymdd = xnc * phi.sin();

        // Euler step
        xm += h * xmd;
        ym += h * ymd;
        xmd += h * xmdd;
        ymd += h * ymdd;
        phi += h * phid;
        t += h;

        // Second derivative evaluation
        let slope = (phibd - phiad) / tf;
        let bint = phibd - slope * tf;
        let mut phid = slope * t + bint;
        if phid > phidmax {
            phid = phidmax;
        }
        if phid < -phidmax {
            phid = -phidmax;
        }
        let xmdd = xnc * phi.cos();
        let ymdd = xnc * phi.sin();

        // RK2 averaging
        xm = (xmold + xm) / 2.0 + 0.5 * h * xmd;
        ym = (ymold + ym) / 2.0 + 0.5 * h * ymd;
        xmd = (xmdold + xmd) / 2.0 + 0.5 * h * xmdd;
        ymd = (ymdold + ymd) / 2.0 + 0.5 * h * ymdd;
        phi = (phiold + phi) / 2.0 + 0.5 * h * phid;

        s += h;
        if s >= ts - 0.0001 {
            s = 0.0;

            let (x1, y1) = predict1(t, xm, ym, xmd, ymd, phi, tf, phiad, phibd, xnc, h, phidmax, accerr);
            let delx = xt - x1;
            let dely = yt - y1;

            let phibd2 = phibd + 0.001;
            let (x2, y2) = predict1(t, xm, ym, xmd, ymd, phi, tf, phiad, phibd2, xnc, h, phidmax, accerr);
            let dxdpb = (x2 - x1) / (phibd2 - phibd);
            let dydpb = (y2 - y1) / (phibd2 - phibd);

            let phiad3 = phiad + 0.001;
            let (x3, y3) = predict1(t, xm, ym, xmd, ymd, phi, tf, phiad3, phibd, xnc, h, phidmax, accerr);
            let dxdpa = (x3 - x1) / (phiad3 - phiad);
            let dydpa = (y3 - y1) / (phiad3 - phiad);

            let (delphiad, delphibd) = if (dxdpb * dydpa - dxdpa * dydpb) == 0.0 {
                (0.0, 0.0)
            } else {
                let delphiad = (dxdpb * dely - delx * dydpb) / (dxdpb * dydpa - dxdpa * dydpb);
                let delphibd = (delx * dydpa - dxdpa * dely) / (dxdpb * dydpa - dxdpa * dydpb);
                (delphiad, delphibd)
            };

            phiad += gain * delphiad;
            phibd += gain * delphibd;

            let phiddeg = phid * 57.3;

            array_t.push(t);
            array_xm.push(xm);
            array_ym.push(ym);
            array_xt.push(xt);
            array_yt.push(yt);
            array_phiddeg.push(phiddeg);
        }
    }

    Results {
        time: array_t,
        xm: array_xm,
        ym: array_ym,
        xt: array_xt,
        yt: array_yt,
        phiddeg: array_phiddeg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c32l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xm.clone(),
        results.ym.clone(),
        results.xt.clone(),
        results.yt.clone(),
        results.phiddeg.clone(),
    ])?;

    let plot_file = format!("{}/c32l4_traj.png", output_dir);
    let config = PlotConfig::new("Rolling Airframe Trajectory (Closed Loop)")
        .with_labels("X (Ft)", "Y (Ft)");

    let series = vec![
        Series::new(results.xm.clone(), results.ym.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
        Series::new(results.xt.clone(), results.yt.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Target"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C32L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c32l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
