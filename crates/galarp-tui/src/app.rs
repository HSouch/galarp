use std::collections::HashMap;
use std::path::{Path, PathBuf};

use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::widgets::ListState;

use crate::configedit::{ConfigDoc, Segment};
use crate::input::LineInput;
use crate::tabs::configs::{self, Check};
use crate::tabs::editor::Editor;


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

/// What mode keystrokes go to.
/// If you're in Normal mode, `q` will quit.
pub enum Mode {
    Normal,
    Help,
    EditField {path: Vec<Segment>, input: LineInput},
}


pub struct App {
    pub quit: bool,
    pub tab: Tab,
    pub mode: Mode,

    pub message: Option<(String, bool)>,

    pub cwd: PathBuf,
    pub files: Vec<PathBuf>,
    pub list: ListState,

    checks: HashMap<PathBuf, Check>,
    pub editor: Option<Editor>,
}


impl App {
    pub fn new(cwd: PathBuf) -> Self {
        let files = configs::find(&cwd);
        let list = ListState::default().with_selected((!files.is_empty()).then_some(0));
        
        Self {
            quit: false,
            tab: Tab::Configs,
            mode: Mode::Normal,
            message: None,
            cwd,
            files,
            list,
            checks: HashMap::new(),
            editor: None,
        }
    }


    /// Turn list index back into a path
    pub fn selected(&self) -> Option<&PathBuf> {
        self.list.selected().and_then(|i| self.files.get(i))
    }


    pub fn selected_check(&mut self) -> Option<&Check> {
        let path = self.selected()?.clone();
        Some(self.checks.entry(path.clone()).or_insert_with(|| configs::check(&path)))
    }


    /// Shows path relative to starting point.
    /// i.e. don't show full path
    pub fn display_path<'a>(&self, p: &'a Path) -> std::borrow::Cow<'a, str> {
        p.strip_prefix(&self.cwd).unwrap_or(p).to_string_lossy()
    }

    fn error(&mut self, msg: impl Into<String>) {
        self.message = Some((msg.into(), true));
    }

    pub fn on_key(&mut self, key: KeyEvent) {
        self.message = None;       // Erase status message on next key press.

        match std::mem::replace(&mut self.mode, Mode::Normal) {
            Mode::Normal => self.on_normal_key(key),
            Mode::Help => { self.mode = Mode::Normal },
            Mode::EditField { path, mut input } => match key.code {
                KeyCode::Esc => {},
                KeyCode::Enter => self.commit_field(&path, &input.text),
                _ => {
                    input.on_key(key);
                    self.mode = Mode::EditField { path, input };
                }
            },
        }
    }

    fn on_normal_key(&mut self, key: KeyEvent) {
        match key.code {
            // q sets mode to quit (will might want to ask for saving later)
            KeyCode::Char('q') => self.quit = true,
            // ? opens a help popup
            KeyCode::Char('?') => self.mode = Mode::Help,
            // Tab and backtab advance tabs around
            KeyCode::Tab => self.tab = Tab::ALL[(self.tab.index() + 1) % Tab::ALL.len()],
            KeyCode::BackTab => self.tab = Tab::ALL[(self.tab.index() + Tab::ALL.len() - 1) % Tab::ALL.len()],
            _ => match self.tab {
                Tab::Configs => self.on_configs_key(key),
                Tab::Edit => self.on_edit_key(key),
                _ => {}
            },
        }
    }    

    fn on_edit_key(&mut self, key: KeyEvent) {
        let Some(ed) = self.editor.as_mut() else { return };

        match key.code {
            KeyCode::Down | KeyCode::Char('j') => ed.move_by(1),
            KeyCode::Up | KeyCode::Char('k') => ed.move_by(-1),
            KeyCode::PageDown => ed.move_by(10),
            KeyCode::PageUp => ed.move_by(-10),
            
            KeyCode::Enter | KeyCode::Char(' ') => {
                let Some(path) = ed.selected_path().map(<[Segment]>::to_vec) else { return };
                match ed.doc.toggle(&path) {
                    Ok(true) => ed.revalidate(),
                    Ok(false) => {
                        let current = ed.doc.value(&path).unwrap_or_default();
                        self.mode = Mode::EditField { path, input: LineInput::new(current)}
                    },
                    Err(e) => self.error(e),
                }
            }

            KeyCode::Delete | KeyCode::Char('d') => {
                let Some(path) = ed.selected_path().map(<[Segment]>::to_vec) else { return };
                self.commit_field(&path, "");
            },
            
            _ => {},
        }
    }


    fn commit_field(&mut self, path: &[Segment], text: &str) {
        let Some(ed) = self.editor.as_mut() else { return };
        match ed.doc.set(path, text) {
            Ok(()) => ed.revalidate(),
            Err(e) => self.error(e),
        }
    }


    fn on_configs_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Down | KeyCode::Char('j') => self.list.select_next(),
            KeyCode::Up | KeyCode::Char('k') => self.list.select_previous(),
            KeyCode::Home | KeyCode::Char('g') => self.list.select_first(),
            KeyCode::End | KeyCode::Char('G') => self.list.select_last(),
            KeyCode::Char('r') => self.rescan(),
            
            KeyCode::Enter | KeyCode::Char('l') | KeyCode::Right => {
                let Some(path) = self.selected().cloned() else { return };
                self.open(&path);
            }

            _ => {},
        }
    }

    /// Loads `path` into the form and switches to the Edit tab.
    fn open(&mut self, path: &Path) {
        match ConfigDoc::load(path) {
            Ok(doc) => {
                self.editor = Some(Editor::new(doc));
                self.tab = Tab::Edit;
            }
            Err(e) => self.error(e),
        }
    }
    

    fn rescan(&mut self) {
        self.files = configs::find(&self.cwd);
        self.checks.clear();
        self.list.select((!self.files.is_empty()).then_some(0));
    }
}
