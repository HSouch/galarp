//! A run to be viewed:
//!     Data, orbit preview, which step is shown, and zoom

use ratatui::crossterm::event::{KeyCode, KeyEvent};

use crate::preview::Preview;
use crate::sim_runner::RunData;


pub struct Viewer {
    pub data: RunData,
    pub preview: Option<Preview>,
    
    /// Log step shown.
    pub cursor: Option<usize>,
    
    // Half-width of the cloud plots (for zooming). If none it is interpreted from the disk size.
    pub view_kpc: Option<f64>,
}


impl Viewer {

    pub fn new(data: RunData, preview: Option<Preview>) -> Self {
        Self { data, preview, cursor: None, view_kpc: None }
    }

    /// Half-width of the cloud plots. This sets the zoom.
    pub fn half_width_kpc(&self) -> f64 {
        self.view_kpc.unwrap_or_else(|| {
            let cfg = &self.data.cfg;
            let region = cfg.diagnostics.disk_region_radius_kpc;
            cfg.galaxy[0].disk.as_ref().map_or(30.0, |d| 1.25 * d.r_max_kpc.max(region))
        })
    }

    pub fn shown_step(&self) -> Option<usize> {
        let n = self.data.steps.len();
        (n > 0).then(|| self.cursor.unwrap_or(n - 1).min(n-1))
    }

    pub fn on_key(&mut self, key: KeyEvent) -> bool {
        let n = self.data.steps.len();

        let at = self.shown_step().unwrap_or(0) as isize;

        // Move to function clamps the move to the bounds set by the steps
        let move_to = |i: isize| Some(i.clamp(0, n.saturating_sub(1) as isize) as usize);

        match key.code {
            KeyCode::Left => self.cursor = move_to(at - 1),
            KeyCode::Right => self.cursor = move_to(at + 1),
            KeyCode::PageUp => self.cursor = move_to(at - 10),
            KeyCode::PageDown => self.cursor = move_to(at + 10),
            KeyCode::Home => self.cursor = Some(0),
            KeyCode::End => self.cursor = None,

            KeyCode::Char('+') | KeyCode::Char('=') => self.view_kpc = Some(self.half_width_kpc() / 1.25),
            KeyCode::Char('-') => self.view_kpc = Some(self.half_width_kpc() * 1.25),
            KeyCode::Char('0') => self.view_kpc = None,
            _ => return false,
        }
        true
    }
}
