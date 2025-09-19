use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct DataXY {
    pub x: Vec<f64>,
    pub y: Vec<f32>,
}

#[derive(Clone, Copy, Debug)]
pub struct FromTo {
    pub from: f64,
    pub to: f64,
}

#[derive(Clone, Debug)]
pub struct Peak {
    pub from: f64,
    pub to: f64,
    pub rt: f64,
    pub integral: f64,
    pub intensity: f64,
    pub ratio: f64,
    pub np: i32,
    pub noise: f64,
}

impl Default for Peak {
    fn default() -> Self {
        Self {
            from: 0.0,
            to: 0.0,
            rt: 0.0,
            integral: 0.0,
            intensity: 0.0,
            ratio: 0.0,
            np: 0,
            noise: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Roi {
    pub rt: f64,
    pub window: f64,
}

impl Roi {
    pub fn new(rt: f64, window: f64) -> Self {
        Self { rt, window }
    }
}

#[derive(Clone, Debug)]
pub struct ChromRoi {
    pub id: String,
    pub idx: usize,
    pub rt: f64,
    pub window: f64,
}

impl ChromRoi {
    pub fn new(id: impl Into<String>, idx: usize, rt: f64, window: f64) -> Self {
        Self {
            id: id.into(),
            idx,
            rt,
            window,
        }
    }
}

#[derive(Clone, Debug)]
pub struct EicRoi {
    pub id: String,
    pub rt: f64,
    pub mz: f64,
    pub window: f64,
}

impl EicRoi {
    pub fn new(id: impl Into<String>, rt: f64, mz: f64, window: f64) -> Self {
        Self {
            id: id.into(),
            rt,
            mz,
            window,
        }
    }
}
