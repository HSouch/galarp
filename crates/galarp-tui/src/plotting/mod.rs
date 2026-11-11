
use ratatui::layout::{Constraint, Rect};
use ratatui::style::{Color, Style, Stylize};
use ratatui::symbols::Marker;
use ratatui::text::Span;
use ratatui::widgets::{Axis, Block, Chart, Dataset, GraphType, LegendPosition};
use ratatui::widgets::canvas::{Canvas, Circle, Line as CanvasLine};
use ratatui::Frame;


pub mod orbits;
pub mod timeseries;

