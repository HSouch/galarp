
use galarp_io::preview;
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


    
}