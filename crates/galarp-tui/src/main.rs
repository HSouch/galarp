

use ratatui::crossterm::event::{ self, Event, KeyEventKind };
use ratatui::DefaultTerminal;


mod app;
mod editor;
mod form;
mod ui;
pub use app::{App, Tab};

// Tabs
mod configs;

use std::time::Duration;




fn main() -> std::io::Result<()> {
    let mut app = app::App::new(std::env::current_dir()?);

    let mut terminal = ratatui::init();
    let res = run(&mut terminal, &mut app);

    ratatui::restore();
    res
}


fn run(terminal: &mut DefaultTerminal, app: &mut app::App ) -> std::io::Result<()> {
    while !app.quit {
        terminal.draw(|frame| ui::draw(frame, app))?;

        if event::poll(Duration::from_millis(50))? {
            if let Event::Key(key) = event::read()? {
                if key.kind == KeyEventKind::Press {
                    app.on_key(key);
                }
            }
        }
    }
    Ok(())
}