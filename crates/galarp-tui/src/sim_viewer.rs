
use ratatui::crossterm::event::{KeyCode, KeyEvent};

use crate::preview::Preview;
use crate::sim_runner::RunData;


pub struct Viewer {
    pub data: RunData,
    pub preview: Option<Preview>,

    pub cursor: Option<usize>,
}


impl Viewer {

    pub fn new(data: RunData) -> Self {
        let preview = Preview::new(&data.cfg).ok();

        Self { data, preview, cursor: None }
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
            _ => return false,
        }

        true
    }

    

}