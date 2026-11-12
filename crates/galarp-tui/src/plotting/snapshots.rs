use super::*;

use ratatui::{text::Line, widgets::canvas::Rectangle};
use ratatui::widgets::canvas::{Points, Line as CanvasLine};

pub const SHADOW_COLORS: [Color; 5] = [
    Color::Blue,
    Color::Cyan,
    Color::Green,
    Color::Yellow,
    Color::Red
];


/// Indication for the `in_disk` diagnostic`
pub enum Region {
    Circle(f64),
    Slab(f64, f64)
}


pub fn cloud_canvas(
    frame: &mut Frame,
    area: Rect,
    title: Line,
    clouds: &[(f64, f64, u8)],
    half: f64,
    region: Region,
    wind: Option<(f64, f64)>,
) {
    let mut bins: [Vec<(f64, f64)>; 5] = Default::default();
    for &(x, y, s) in clouds {
        bins[(s as usize * 5 / 256).min(4)].push((x, y));
    }   

    let (hw, hh) = square_half_sizes(area, half);

    let canvas = Canvas::default()
        .block(Block::bordered().title(title))
        .marker(Marker::Braille)
        .x_bounds([-hw, hw])
        .y_bounds([-hh, hh])
        .paint(|ctx| {
            match region {
                Region::Circle(r) => ctx.draw(&Circle {x: 0.0, y: 0.0, radius: r, color: Color::DarkGray}),
                Region::Slab(r, h) => {
                    ctx.draw(&Rectangle {x: -r, y: -h, width:2.0 * r, height: 2.0 * h, color: Color::DarkGray})
                }
            }
            ctx.layer();
            for (pts, &color) in bins.iter().zip(&SHADOW_COLORS) {
                ctx.draw(&Points { coords: pts, color});
            }
            if let Some((wx, wy)) = wind {
                let (len, x0, y0) = (0.25 * hw.min(hh), -0.77 * hw, 0.75 * hh);
                ctx.layer();
                ctx.draw(&CanvasLine { x1: x0, y1: y0, x2: x0 + len * wx, y2: y0 + len * wy, color: Color::White});

                ctx.print(x0 + len * wx, y0 + len * wy, "*".white());
                ctx.print(x0 -0.1 * hw, y0 + 0.12 * hh, "wind".white().dim());
            }
        });

    frame.render_widget(canvas, area);
}