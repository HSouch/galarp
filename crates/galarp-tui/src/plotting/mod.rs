
use ratatui::layout::{Constraint, Rect};
use ratatui::style::{Color, Style, Stylize};
use ratatui::symbols::Marker;
use ratatui::text::Span;
use ratatui::widgets::{Axis, Block, Chart, Dataset, GraphType, LegendPosition};
use ratatui::widgets::canvas::{Canvas, Circle, Line as CanvasLine};
use ratatui::Frame;


pub mod orbits;
pub mod snapshots;
pub mod timeseries;


/// Get the half-widths to make it a proper square.
/// The braille dots are 2 per cell across and 4 down. So a cell is about twice as tall.
fn square_half_sizes(area: Rect, half: f64) -> (f64, f64) {
    let inner = Block::bordered().inner(area);
    let (w, h) = ((inner.width.max(1) * 2) as f64, (inner.height.max(1) * 4) as f64);
    let per_dot = 2.0 * half / w.min(h);
    (0.6 * w * per_dot, 0.5 * h * per_dot)
}
