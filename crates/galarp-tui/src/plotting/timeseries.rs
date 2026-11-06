
// use ratatui::layout::Rect;
// use ratatui::style::{Color, Style};
// use ratatui::symbols::Marker;
// use ratatui::text::Span;
// use ratatui::widgets::{Axis, Block, Chart, Dataset, GraphType};
// use ratatui::Frame;

use super::*;

#[derive(Clone, Copy)]
pub enum YAxis {
    /// Fixed range
    Fixed(f64, f64),
    /// Fitted to data
    Auto,
}

pub fn time_xy_chart(frame: &mut Frame, area: Rect, title: &str, data: &[(f64, f64)], t_end: f64,
    y: YAxis, cursor: Option<f64>) {
    
    let finite: Vec<(f64, f64)> = data.iter().copied()
        .filter(|p| p.1.is_finite()).collect();

    let (mut lo, mut hi) = match y {
        YAxis::Fixed(lo, hi) => (lo, hi),
        YAxis::Auto => finite.iter()
            .fold(
                (f64::INFINITY, f64::NEG_INFINITY),
                |(a, b), p| (a.min(p.1), b.max(p.1))
            ),
    };

    /// If auto, log the y-axis if it spans two orders of magnitude
    let log = matches!(y, YAxis::Auto) && lo > 0.0 && hi / lo > 100.0;
    let points: Vec<(f64, f64)> = if log { finite.iter().map(|&(t, v)| (t, v.log10())).collect() } else {finite };


    if log {
        (lo, hi) = (lo.log10(), hi.log10());
    }
    if !(lo.is_finite() && hi.is_finite()) {
        (lo, hi) = (0.0, 1.0);
    }

    if hi - lo < 1e-12 * hi.abs().max(1e-300) {
        let pad = if hi == 0.0 { 1.0 } else { 0.1 * hi.abs() };
        (lo, hi) = (lo - pad, hi + pad);
    }

    
    // let label = |v: f64| if log { fmt(10f64.powf(v)) } else { fmt(v) };

    let line = |color: Color| 
        Dataset::default().marker(Marker::Braille).graph_type(GraphType::Line)
            .style(Style::new().fg(color));

    let mut datasets = vec![
        line(Color::Cyan).data(&points)
    ];

    let chart = Chart::new(datasets)
        .block(Block::bordered().title(title))
        .x_axis(Axis::default().bounds([0.0, t_end]).style(Style::new().dim()))
        .y_axis(Axis::default().bounds([lo, hi]).style(Style::new().dim()));

    frame.render_widget(chart, area);

}