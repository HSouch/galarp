
use super::*;

#[derive(Clone, Copy)]
pub enum YAxis {
    /// Fixed range
    Fixed(f64, f64),
    /// Fitted to data
    Auto,
}

/// One line series for charting
pub struct Series<'a> {
    pub name: &'a str,
    pub color: Color,
    pub data: Vec<(f64, f64)>,
}


pub fn time_xy_chart(frame: &mut Frame, area: Rect, title: &str, data: &[(f64, f64)], t_end: f64,
    y: YAxis, cursor: Option<f64>) {
        /// Line chart of (t, y)
        
        let series = [Series {name: "", color: Color::Cyan, data: data.to_vec()}];
        chart(frame, area, title, &series, t_end, y, cursor);
    }


pub fn chart(frame: &mut Frame, area: Rect, title: &str, series: &[Series], t_end: f64, y: YAxis, cursor: Option<f64>) {
    
    let finite: Vec<Vec<(f64, f64)>> = 
        series.iter().map(|s| s.data.iter().copied().filter(
            |p| p.1.is_finite()).collect()).collect();
    
    
    let (mut lo, mut hi) = match y {
        YAxis::Fixed(lo, hi) => (lo, hi),
        YAxis::Auto => finite.iter().flatten()
            .fold(
                (f64::INFINITY, f64::NEG_INFINITY),
                |(a, b), p| (a.min(p.1), b.max(p.1))
            ),
    };

    // If auto, log the y-axis if it spans two orders of magnitude
    let log = matches!(y, YAxis::Auto) && lo > 0.0 && hi / lo > 100.0;
    let points: Vec<Vec<(f64, f64)>> = if log { finite.iter().map(
        |d| d.iter().map(
            |&(t, v)| (t, v.log10())).collect()).collect()
        } else {finite };

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

    let cursor_line: Vec<(f64, f64)> = cursor.map(|t| vec![(t, lo), (t, hi)]).unwrap_or_default();

    let line = |color: Color|
        Dataset::default()
        .marker(Marker::Braille)
        .graph_type(GraphType::Line)
        .style(Style::new().fg(color));

    let mut datasets: Vec<Dataset> = series
        .iter()
        .zip(&points)
        .map(|(s, d)| {
            let ds = line(s.color).data(d);
            if series.len() > 1 { ds.name(s.name) } else { ds }
        }).collect();
    
    if cursor.is_some() {
        datasets.push(line(Color::Yellow).data(&cursor_line));
    }

    let chart = Chart::new(datasets)
        .block(Block::bordered().title(title))
        .x_axis(Axis::default().bounds([0.0, t_end]).style(Style::new().dim()))
        .y_axis(Axis::default().bounds([lo, hi]).style(Style::new().dim()));

    frame.render_widget(chart, area);

}