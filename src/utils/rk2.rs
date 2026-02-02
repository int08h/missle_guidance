//! Runge-Kutta 2nd order integration helper
//!
//! Provides generic RK2 integration utilities used throughout the simulations

/// Generic state type for numerical integration
pub trait State: Clone {
    /// Add two states element-wise
    fn add(&self, other: &Self) -> Self;

    /// Multiply state by a scalar
    fn scale(&self, factor: f64) -> Self;
}

/// Implement State for arrays of different sizes
impl State for [f64; 2] {
    fn add(&self, other: &Self) -> Self {
        [self[0] + other[0], self[1] + other[1]]
    }

    fn scale(&self, factor: f64) -> Self {
        [self[0] * factor, self[1] * factor]
    }
}

impl State for [f64; 4] {
    fn add(&self, other: &Self) -> Self {
        [
            self[0] + other[0],
            self[1] + other[1],
            self[2] + other[2],
            self[3] + other[3],
        ]
    }

    fn scale(&self, factor: f64) -> Self {
        [
            self[0] * factor,
            self[1] * factor,
            self[2] * factor,
            self[3] * factor,
        ]
    }
}

impl State for [f64; 6] {
    fn add(&self, other: &Self) -> Self {
        [
            self[0] + other[0],
            self[1] + other[1],
            self[2] + other[2],
            self[3] + other[3],
            self[4] + other[4],
            self[5] + other[5],
        ]
    }

    fn scale(&self, factor: f64) -> Self {
        [
            self[0] * factor,
            self[1] * factor,
            self[2] * factor,
            self[3] * factor,
            self[4] * factor,
            self[5] * factor,
        ]
    }
}

impl State for Vec<f64> {
    fn add(&self, other: &Self) -> Self {
        self.iter().zip(other.iter()).map(|(a, b)| a + b).collect()
    }

    fn scale(&self, factor: f64) -> Self {
        self.iter().map(|x| x * factor).collect()
    }
}

/// RK2 integrator matching the MATLAB pattern used in the simulations
///
/// This implements the specific RK2 variant used throughout the book:
/// 1. Euler step forward
/// 2. Evaluate derivatives at new position
/// 3. Average (old + new)/2 + 0.5*h*derivative
pub struct RK2Integrator<S: State, F: Fn(&S, f64) -> S> {
    pub state: S,
    pub time: f64,
    pub step: f64,
    derivative_fn: F,
}

impl<S: State, F: Fn(&S, f64) -> S> RK2Integrator<S, F> {
    pub fn new(initial_state: S, initial_time: f64, step_size: f64, derivative_fn: F) -> Self {
        Self {
            state: initial_state,
            time: initial_time,
            step: step_size,
            derivative_fn,
        }
    }

    /// Take one RK2 step
    pub fn step(&mut self) {
        let state_old = self.state.clone();
        let h = self.step;

        // First derivative evaluation
        let k1 = (self.derivative_fn)(&self.state, self.time);

        // Euler step
        self.state = self.state.add(&k1.scale(h));
        self.time += h;

        // Second derivative evaluation
        let k2 = (self.derivative_fn)(&self.state, self.time);

        // RK2 correction: (old + euler_result)/2 + 0.5*h*k2
        let sum = state_old.add(&self.state);
        let avg = sum.scale(0.5);
        self.state = avg.add(&k2.scale(0.5 * h));
    }

    /// Integrate to a target time
    pub fn integrate_to(&mut self, target_time: f64) {
        while self.time < target_time - 0.00001 {
            self.step();
        }
    }
}

/// Simple 2D vector operations
#[derive(Debug, Clone, Copy, Default)]
pub struct Vec2 {
    pub x: f64,
    pub y: f64,
}

impl Vec2 {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();
        if mag > 0.0 {
            Self { x: self.x / mag, y: self.y / mag }
        } else {
            Self { x: 0.0, y: 0.0 }
        }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y
    }
}

impl std::ops::Add for Vec2 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self { x: self.x + other.x, y: self.y + other.y }
    }
}

impl std::ops::Sub for Vec2 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self { x: self.x - other.x, y: self.y - other.y }
    }
}

impl std::ops::Mul<f64> for Vec2 {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self { x: self.x * scalar, y: self.y * scalar }
    }
}

/// Simple 3D vector operations
#[derive(Debug, Clone, Copy, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();
        if mag > 0.0 {
            Self { x: self.x / mag, y: self.y / mag, z: self.z / mag }
        } else {
            Self { x: 0.0, y: 0.0, z: 0.0 }
        }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}

impl std::ops::Add for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self { x: self.x + other.x, y: self.y + other.y, z: self.z + other.z }
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self { x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
    }
}

impl std::ops::Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        Self { x: self.x * scalar, y: self.y * scalar, z: self.z * scalar }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec2_operations() {
        let v1 = Vec2::new(3.0, 4.0);
        assert!((v1.magnitude() - 5.0).abs() < 1e-10);

        let v2 = Vec2::new(1.0, 0.0);
        assert!((v1.dot(&v2) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_vec3_cross() {
        let i = Vec3::new(1.0, 0.0, 0.0);
        let j = Vec3::new(0.0, 1.0, 0.0);
        let k = i.cross(&j);
        assert!((k.z - 1.0).abs() < 1e-10);
    }
}
