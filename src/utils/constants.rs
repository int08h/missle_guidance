//! Physical constants used throughout the simulations

use std::f64::consts::PI;

/// Earth's radius in feet
pub const EARTH_RADIUS_FT: f64 = 2.0926e7;

/// Earth's radius in kilometers
pub const EARTH_RADIUS_KM: f64 = 6380.0;

/// Earth's gravitational parameter (GM) in ft^3/sec^2
pub const GM_FT: f64 = 1.4077e16;

/// Earth's gravitational parameter (GM) in km^3/sec^2
pub const GM_KM: f64 = 398923.0;

/// Gravitational acceleration at sea level (ft/s^2)
pub const G_ACCEL: f64 = 32.2;

/// Pi constant for precision matching with MATLAB
pub const PI_VAL: f64 = 3.1415926535898;

/// Degrees per radian
pub const DEG_PER_RAD: f64 = 57.29577951308232;

/// Radians per degree
pub const RAD_PER_DEG: f64 = 1.0 / DEG_PER_RAD;

/// Half PI
pub const HALF_PI: f64 = PI / 2.0;

/// Feet per kilometer
pub const FT_PER_KM: f64 = 3280.0;

/// Convert degrees to radians
#[inline]
pub fn deg_to_rad(deg: f64) -> f64 {
    deg / DEG_PER_RAD
}

/// Convert radians to degrees
#[inline]
pub fn rad_to_deg(rad: f64) -> f64 {
    rad * DEG_PER_RAD
}
