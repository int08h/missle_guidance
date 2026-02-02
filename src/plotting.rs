//! Plotting utilities for simulation visualization
//!
//! Uses the plotters crate to generate plots similar to MATLAB output

use plotters::prelude::*;

/// Plot configuration
pub struct PlotConfig {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub width: u32,
    pub height: u32,
    pub x_range: Option<(f64, f64)>,
    pub y_range: Option<(f64, f64)>,
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self {
            title: String::new(),
            x_label: "X".to_string(),
            y_label: "Y".to_string(),
            width: 800,
            height: 600,
            x_range: None,
            y_range: None,
        }
    }
}

impl PlotConfig {
    pub fn new(title: &str) -> Self {
        Self {
            title: title.to_string(),
            ..Default::default()
        }
    }

    pub fn with_labels(mut self, x_label: &str, y_label: &str) -> Self {
        self.x_label = x_label.to_string();
        self.y_label = y_label.to_string();
        self
    }

    pub fn with_x_range(mut self, min: f64, max: f64) -> Self {
        self.x_range = Some((min, max));
        self
    }

    pub fn with_y_range(mut self, min: f64, max: f64) -> Self {
        self.y_range = Some((min, max));
        self
    }

    #[allow(dead_code)]
    pub fn with_size(mut self, width: u32, height: u32) -> Self {
        self.width = width;
        self.height = height;
        self
    }
}

/// Series data for plotting
pub struct Series {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub label: Option<String>,
    pub color: RGBColor,
}

impl Series {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        Self {
            x,
            y,
            label: None,
            color: BLUE,
        }
    }

    pub fn with_label(mut self, label: &str) -> Self {
        self.label = Some(label.to_string());
        self
    }

    pub fn with_color(mut self, color: RGBColor) -> Self {
        self.color = color;
        self
    }
}

/// Create a line plot and save to PNG
pub fn line_plot(filename: &str, config: &PlotConfig, series: &[Series]) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (config.width, config.height)).into_drawing_area();
    root.fill(&WHITE)?;

    // Calculate ranges
    let (x_min, x_max) = config.x_range.unwrap_or_else(|| {
        let all_x: Vec<f64> = series.iter().flat_map(|s| s.x.iter().copied()).collect();
        let min = all_x.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = all_x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let margin = (max - min) * 0.05;
        (min - margin, max + margin)
    });

    let (y_min, y_max) = config.y_range.unwrap_or_else(|| {
        let all_y: Vec<f64> = series.iter().flat_map(|s| s.y.iter().copied()).collect();
        let min = all_y.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = all_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let margin = (max - min) * 0.05;
        (min - margin, max + margin)
    });

    let mut chart = ChartBuilder::on(&root)
        .caption(&config.title, ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc(&config.x_label)
        .y_desc(&config.y_label)
        .draw()?;

    for s in series {
        let data: Vec<(f64, f64)> = s.x.iter().zip(s.y.iter()).map(|(&x, &y)| (x, y)).collect();

        if let Some(ref label) = s.label {
            chart
                .draw_series(LineSeries::new(data, s.color.stroke_width(2)))?
                .label(label)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], s.color.stroke_width(2)));
        } else {
            chart.draw_series(LineSeries::new(data, s.color.stroke_width(2)))?;
        }
    }

    if series.iter().any(|s| s.label.is_some()) {
        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.8))
            .border_style(BLACK)
            .draw()?;
    }

    root.present()?;
    Ok(())
}

/// Create a scatter plot
pub fn scatter_plot(
    filename: &str,
    config: &PlotConfig,
    x: &[f64],
    y: &[f64],
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (config.width, config.height)).into_drawing_area();
    root.fill(&WHITE)?;

    let (x_min, x_max) = config.x_range.unwrap_or_else(|| {
        let min = x.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let margin = (max - min) * 0.05;
        (min - margin, max + margin)
    });

    let (y_min, y_max) = config.y_range.unwrap_or_else(|| {
        let min = y.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let margin = (max - min) * 0.05;
        (min - margin, max + margin)
    });

    let mut chart = ChartBuilder::on(&root)
        .caption(&config.title, ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_desc(&config.x_label)
        .y_desc(&config.y_label)
        .draw()?;

    let data: Vec<(f64, f64)> = x.iter().zip(y.iter()).map(|(&x, &y)| (x, y)).collect();
    chart.draw_series(
        data.iter().map(|&(x, y)| Circle::new((x, y), 3, BLUE.filled())),
    )?;

    root.present()?;
    Ok(())
}

/// Quick plot with minimal configuration
#[allow(dead_code)]
pub fn quick_plot(
    filename: &str,
    x: &[f64],
    y: &[f64],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let config = PlotConfig::new(title).with_labels(x_label, y_label);
    let series = vec![Series::new(x.to_vec(), y.to_vec())];
    line_plot(filename, &config, &series)
}

/// Multiple series quick plot
#[allow(dead_code)]
pub fn quick_multi_plot(
    filename: &str,
    data: &[(Vec<f64>, Vec<f64>, &str)],
    title: &str,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let colors = [BLUE, RED, GREEN, MAGENTA, CYAN, BLACK];
    let config = PlotConfig::new(title).with_labels(x_label, y_label);

    let series: Vec<Series> = data
        .iter()
        .enumerate()
        .map(|(i, (x, y, label))| {
            Series::new(x.clone(), y.clone())
                .with_label(label)
                .with_color(colors[i % colors.len()])
        })
        .collect();

    line_plot(filename, &config, &series)
}
