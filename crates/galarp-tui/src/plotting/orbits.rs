
use ratatui::{widgets::canvas::{Canvas, Circle, Line as CanvasLine}};
use galarp::utils::{dot};

use super::*;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Projection {
    XY,
    XZ,
    YZ,
}


impl Projection {
    pub fn next(self) -> Self {
        match self {
            Projection::XY => Projection::XZ,
            Projection::XZ => Projection::YZ,
            Projection::YZ => Projection::XY,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Projection::XY => "x-y",
            Projection::XZ => "x-z",
            Projection::YZ => "y-z",
        }
    }

    fn axes(self) -> ([f64; 3], [f64;3]) {
        match self {
            Projection::XY => ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            Projection::XZ => ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            Projection::YZ => ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        }
    }
}


pub fn orbit_canvas(frame: &mut Frame, area: Rect, path: &[[f64; 3]], r_vir: Option<f64>,
    proj: Projection, current: Option<usize>) {
        
        let (e1, e2) = proj.axes();

        let pts: Vec<(f64, f64)> = path.iter().map(|p| (dot(p, &e1), dot(p, &e2))).collect();
        
        // Get the boundaries for the orbit
        let (mut x0, mut x1, mut y0, mut y1) = pts.iter()
            .fold((0.0f64, 0.0f64, 0.0f64, 0.0f64),
            |(a, b, c, d), p| (a.min(p.0), b.max(p.0), c.min(p.1), d.max(p.1)));

        if let Some(r) = r_vir {
            (x0, x1, y0, y1) = (x0.min(-r), x1.max(r), y0.min(-r), y1.max(r))
        }

        let span = (x1 - x0).max(y1 - y0).max(1.0) * 1.1;
        let (cx, cy) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
        let (hw, hh) = square_half_sizes(area, 0.5 * span);

        let title = format!(" orbit ({}), {} kpc across", proj.name(), 2.0 * hw);

        let canvas = Canvas::default()
            .block(Block::bordered().title(title))
            .marker(Marker::Braille)
            .x_bounds([cx - hw, cx + hw])
            .y_bounds([cy - hh, cy + hh])
            .paint( |ctx| {
                // Draw the beta profile core radius
                if let Some(r) = r_vir {
                    ctx.draw(&Circle { x: 0.0, y: 0.0, radius: r, color: Color::DarkGray});
                }
                // Draw the orbit through the cluster
                for w in pts.windows(2) {
                    ctx.draw(&CanvasLine { x1: w[0].0, y1: w[0].1, x2: w[1].0, y2: w[1].1, color: Color::LightMagenta});
                }

                if let Some(&(x, y)) = current.and_then(|i| pts.get(i)) {
                    ctx.print(x, y, "X");
                }
            });

        frame.render_widget(canvas, area);
    }

