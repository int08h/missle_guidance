//! Chapter 39, Lesson 3: Optimal Guidance with Full Autopilot Model
//!
//! Full autopilot model simulation with optimal guidance.

use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub miss: Vec<f64>,
}

/// Run the C39L3 simulation
pub fn run() -> Results {
    // System parameters
    let alt: f64 = 100000.0;
    let tau: f64 = 0.5;
    let xnt: f64 = 161.0;
    let xlim: f64 = 322.0;
    let apn: i32 = 2;
    let wcr: f64 = 50.0;
    let h: f64 = 0.001;
    let xnp: f64 = 3.0;
    let wact: f64 = 150.0;
    let zact: f64 = 0.7;
    let vc: f64 = 4000.0;
    let vm: f64 = 3000.0;
    let fr: f64 = 3.0;
    let diam: f64 = 1.0;
    let xl: f64 = 20.0;
    let ctw: f64 = 0.0;
    let crw: f64 = 6.0;
    let hw: f64 = 2.0;
    let ctt: f64 = 0.0;
    let crt: f64 = 2.0;
    let ht: f64 = 2.0;
    let xn: f64 = 4.0;
    let xcg: f64 = 10.0;
    let xhl: f64 = 19.5;
    let zeta: f64 = 0.7;
    let a: f64 = 1000.0;
    let gam: f64 = 0.001;

    // Atmosphere
    let rho = if alt <= 30000.0 {
        0.002378 * (-alt / 30000.0).exp()
    } else {
        0.0034 * (-alt / 22000.0).exp()
    };

    let wgt: f64 = 1000.0;
    let _xacc = xcg;
    let swing = 0.5 * hw * (ctw + crw);
    let stail = 0.5 * ht * (ctt + crt);
    let sref = std::f64::consts::PI * diam * diam / 4.0;
    let xlp = fr * diam;
    let splan = (xl - xlp) * diam + 1.33 * xlp * diam / 2.0;
    let xcpn = 2.0 * xlp / 3.0;
    let an = 0.67 * xlp * diam;
    let ab = (xl - xlp) * diam;
    let xcpb = (0.67 * an * xlp + ab * (xlp + 0.5 * (xl - xlp))) / (an + ab);
    let xcpw = xlp + xn + 0.7 * crw - 0.2 * ctw;
    let xmach = vm / a;
    let xiyy = wgt * (3.0 * ((diam / 2.0).powi(2)) + xl * xl) / (12.0 * 32.2);

    let tmp1 = (xcg - xcpw) / diam;
    let tmp2 = (xcg - xhl) / diam;
    let tmp3 = (xcg - xcpb) / diam;
    let tmp4 = (xcg - xcpn) / diam;
    let b = (xmach.powi(2) - 1.0).sqrt();
    let q = 0.5 * rho * vm * vm;
    let xncg = 5.0;
    let p1 = wgt * xncg / (q * sref);
    let y1 = 2.0 + 8.0 * swing / (b * sref) + 8.0 * stail / (b * sref);
    let y2 = 1.5 * splan / sref;
    let y3 = 8.0 * stail / (b * sref);
    let y4 = 2.0 * tmp4 + 8.0 * swing * tmp1 / (b * sref) + 8.0 * stail * tmp2 / (b * sref);
    let y5 = 1.5 * splan * tmp3 / sref;
    let y6 = 8.0 * stail * tmp2 / (b * sref);
    let p2 = y2 - y3 * y5 / y6;
    let p3 = y1 - y3 * y4 / y6;
    let alftr = (-p3 + (p3 * p3 + 4.0 * p2 * p1).sqrt()) / (2.0 * p2);
    let _deltr = -y4 * alftr / y6 - y5 * alftr * alftr / y6;
    let cna = 2.0 + 1.5 * splan * alftr / sref + 8.0 * swing / (b * sref) + 8.0 * stail / (b * sref);
    let cnd = 8.0 * stail / (b * sref);
    let cmap = 2.0 * tmp4 + 1.5 * splan * alftr * tmp3 / sref + 8.0 * swing * tmp1 / (b * sref);
    let cma = cmap + 8.0 * stail * tmp2 / (b * sref);
    let cmd = 8.0 * stail * tmp2 / (b * sref);
    let xma = q * sref * diam * cma / xiyy;
    let xmd = q * sref * diam * cmd / xiyy;
    let za = -32.2 * q * sref * cna / (wgt * vm);
    let zd = -32.2 * q * sref * cnd / (wgt * vm);
    let wz = ((xma * zd - za * xmd) / zd).sqrt();
    let waf = (-xma).sqrt();
    let zaf = 0.5 * waf * za / xma;
    let xk1 = -vm * (xma * zd - xmd * za) / (1845.0 * xma);
    let _xk2 = xk1;
    let ta = xmd / (xma * zd - xmd * za);
    let xk3 = 1845.0 * xk1 / vm;
    let w = (tau * wcr * (1.0 + 2.0 * zaf * waf / wcr) - 1.0) / (2.0 * zeta * tau);
    let w0 = w / (tau * wcr).sqrt();
    let z0 = 0.5 * w0 * (2.0 * zeta / w + tau - waf.powi(2) / (w0 * w0 * wcr));
    let xkc = (-w0.powi(2) / wz.powi(2) - 1.0 + 2.0 * z0 * w0 * ta) / (1.0 - 2.0 * z0 * w0 * ta + w0 * w0 * ta * ta);
    let xka = xk3 / (xk1 * xkc);
    let xk0 = -w * w / (tau * waf * waf);
    let xk = xk0 / (xk1 * (1.0 + xkc));
    let wi = xkc * ta * w0 * w0 / (1.0 + xkc + w0.powi(2) / wz.powi(2));
    let xkr = xk / (xka * wi);
    let xkdc = 1.0 + 1845.0 / (xka * vm);

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf = 0.1;
    while tf <= 10.0 {
        let mut x: f64 = 0.0;
        let mut t: f64 = 0.0;
        let mut _s: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut yd: f64 = 0.0;
        let _xnc: f64 = 0.0;
        let mut xnl: f64 = 0.0;
        let mut deld: f64 = 0.0;
        let mut del: f64 = 0.0;
        let mut e: f64 = 0.0;
        let mut ed: f64 = 0.0;

        while t <= tf - 1e-5 {
            _s += h;
            let xold = x;
            let yold = y;
            let ydold = yd;
            let delold = del;
            let deldold = deld;
            let eold = e;
            let edold = ed;

            // First derivative evaluation
            let tgo = tf - t + 0.0001;
            let rtm = vc * tgo;
            let _xlam = y / rtm;
            let xlamd = (y + yd * tgo) / (vc * tgo * tgo);

            let mut xnc_new = if apn == 0 {
                xnp * vc * xlamd
            } else if apn == 1 {
                let xs = tgo / tau;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs * xs * xs + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let xnlpz = xnl * 32.2;
                let xnew = xnpp * xnlpz * ((-xs).exp() + xs - 1.0) / (xs * xs);
                xnpp * vc * xlamd + 0.5 * xnpp * xnt - xnew
            } else {
                let xs = tgo / tau;
                let temp1 = tgo * tgo * tau * ((-xs).exp() - 1.0 + xs);
                let top = -(tgo.powi(3)) / (tau * wz) + (1.0 + 1.0 / (tau * wz)) * temp1;
                let temp2 = 0.5 * (1.0 - 3.0 / (tau * wz)) + xs * (1.0 + 1.0 / (tau * wz)) - xs * xs;
                let temp3 = -2.0 * xs * (-xs).exp();
                let temp4 = 2.0 * (-xs).exp() / (tau * wz) - 0.5 * (-2.0 * xs).exp() * (1.0 + 1.0 / (tau * wz));
                let bot = gam + tgo.powi(3) / 3.0 + (1.0 + 1.0 / (tau * wz)) * tau.powi(3) * (temp2 + temp3 + temp4);
                let xnpp = top / bot;
                let c1th = xnpp / (tgo * tgo);
                let c2th = xnpp / tgo;
                let c3th = 0.5 * xnpp;
                let c4th = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                let xnlpz = xnl * 32.2;
                let epz = xnlpz + 1.0 / (tau * wz);
                c1th * y + c2th * yd + c3th * xnt + c4th * epz
            };

            xnc_new = xnc_new.clamp(-xlim, xlim);
            let xncg = xnc_new / 32.2;

            let thd = xk3 * (e + ta * ed);
            let delc = xkr * (x + thd);
            let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
            let edd = waf * waf * (del - e - 2.0 * zaf * ed / waf);
            xnl = xk1 * (e - edd / (wz * wz));
            let xd = wi * (thd + xka * (xnl - xncg * xkdc));
            let xnlpz = xnl * 32.2;
            let ydd = xnt - xnlpz;

            // Euler step
            x += h * xd;
            y += h * yd;
            yd += h * ydd;
            del += h * deld;
            deld += h * deldd;
            e += h * ed;
            ed += h * edd;
            t += h;

            // Second derivative evaluation for RK2
            let tgo = tf - t + 0.0001;
            let _rtm = vc * tgo;
            let _xlam = y / (vc * tgo);
            let xlamd = (y + yd * tgo) / (vc * tgo * tgo);

            let mut xnc_new = if apn == 0 {
                xnp * vc * xlamd
            } else if apn == 1 {
                let xs = tgo / tau;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs * xs * xs + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let xnlpz = xnl * 32.2;
                let xnew = xnpp * xnlpz * ((-xs).exp() + xs - 1.0) / (xs * xs);
                xnpp * vc * xlamd + 0.5 * xnpp * xnt - xnew
            } else {
                let xs = tgo / tau;
                let temp1 = tgo * tgo * tau * ((-xs).exp() - 1.0 + xs);
                let top = -(tgo.powi(3)) / (tau * wz) + (1.0 + 1.0 / (tau * wz)) * temp1;
                let temp2 = 0.5 * (1.0 - 3.0 / (tau * wz)) + xs * (1.0 + 1.0 / (tau * wz)) - xs * xs;
                let temp3 = -2.0 * xs * (-xs).exp();
                let temp4 = 2.0 * (-xs).exp() / (tau * wz) - 0.5 * (-2.0 * xs).exp() * (1.0 + 1.0 / (tau * wz));
                let bot = gam + tgo.powi(3) / 3.0 + (1.0 + 1.0 / (tau * wz)) * tau.powi(3) * (temp2 + temp3 + temp4);
                let xnpp = top / bot;
                let c1th = xnpp / (tgo * tgo);
                let c2th = xnpp / tgo;
                let c3th = 0.5 * xnpp;
                let c4th = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                let xnlpz = xnl * 32.2;
                let epz = xnlpz + 1.0 / (tau * wz);
                c1th * y + c2th * yd + c3th * xnt + c4th * epz
            };

            xnc_new = xnc_new.clamp(-xlim, xlim);
            let xncg = xnc_new / 32.2;

            let thd = xk3 * (e + ta * ed);
            let delc = xkr * (x + thd);
            let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
            let edd = waf * waf * (del - e - 2.0 * zaf * ed / waf);
            xnl = xk1 * (e - edd / (wz * wz));
            let xd = wi * (thd + xka * (xnl - xncg * xkdc));
            let xnlpz = xnl * 32.2;
            let ydd = xnt - xnlpz;

            // RK2 averaging
            x = 0.5 * (xold + x + h * xd);
            y = 0.5 * (yold + y + h * yd);
            yd = 0.5 * (ydold + yd + h * ydd);
            del = 0.5 * (delold + del + h * deld);
            deld = 0.5 * (deldold + deld + h * deldd);
            e = 0.5 * (eold + e + h * ed);
            ed = 0.5 * (edold + ed + h * edd);
        }

        array_tf.push(tf);
        array_y.push(y);
        tf += 0.1;
    }

    Results {
        tf: array_tf,
        miss: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c39l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.miss.clone(),
    ])?;

    println!("C39L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c39l3_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
