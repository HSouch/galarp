
use ratatui::crossterm::event::{KeyCode, KeyEvent};


#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Configs,
    Edit,
    Run,
    Results,
    Wizard,
}


impl Tab {
    pub const ALL: [Tab; 5] = [
        Tab::Configs,
        Tab::Edit,
        Tab::Run,
        Tab::Results,
        Tab::Wizard,
    ];

    pub fn title(self) -> &'static str {
        match self {
            Tab::Configs => "Configs",
            Tab::Edit => "Edit",
            Tab::Run => "Run",
            Tab::Results => "Results",
            Tab::Wizard => "Wizard",
        }
    }

    fn index(self) -> usize {
        Self::ALL.iter().position(|&t| t == self).unwrap()
    }
}


pub struct App {
    pub quit: bool,
    pub tab: Tab,
}


impl App {
    pub fn new() -> Self {
        Self {
            quit: false,
            tab: Tab::Configs,
        }
    }

    pub fn ok_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Char('q') => self.quit = true,
            KeyCode::Tab => self.tab = Tab::ALL[(self.tab.index() + 1) % Tab::ALL.len()],
            KeyCode::BackTab => self.tab = Tab::ALL[(self.tab.index() + Tab::ALL.len() - 1) % Tab::ALL.len()],
            _ => {},
        }
    }
}