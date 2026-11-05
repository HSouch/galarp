
use std::path::{Path, PathBuf};

use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::widgets::ListState;

use crate::tabs::configs;


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

    pub cwd: PathBuf,
    pub files: Vec<PathBuf>,
    pub list: ListState,
}


impl App {
    pub fn new(cwd: PathBuf) -> Self {
        let files = configs::find(&cwd);
        let list = ListState::default().with_selected((!files.is_empty()).then_some(0));
        
        Self {
            quit: false,
            tab: Tab::Configs,
            cwd,
            files,
            list
        }
    }


    /// Turn list index back into a path
    pub fn selected(&self) -> Option<&PathBuf> {
        self.list.selected().and_then(|i| self.files.get(i))
    }


    /// Shows path relative to starting point.
    /// i.e. don't show full path
    pub fn display_path<'a>(&self, p: &'a Path) -> std::borrow::Cow<'a, str> {
        p.strip_prefix(&self.cwd).unwrap_or(p).to_string_lossy()
    }


    pub fn on_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Char('q') => self.quit = true,
            KeyCode::Tab => self.tab = Tab::ALL[(self.tab.index() + 1) % Tab::ALL.len()],
            KeyCode::BackTab => self.tab = Tab::ALL[(self.tab.index() + Tab::ALL.len() - 1) % Tab::ALL.len()],
            _ => match self.tab {
                Tab::Configs => self.on_configs_key(key),
                _ => {}
            },
        }
    }


    fn on_configs_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Down | KeyCode::Char('j') => self.list.select_next(),
            KeyCode::Up | KeyCode::Char('k') => self.list.select_previous(),
            KeyCode::Home | KeyCode::Char('g') => self.list.select_first(),
            KeyCode::End | KeyCode::Char('G') => self.list.select_last(),
            KeyCode::Char('r') => self.rescan(),
            _ => {},
        }
    }
    

    fn rescan(&mut self) {
        self.files = configs::find(&self.cwd);
        self.list.select((!self.files.is_empty()).then_some(0));
    }
}
