

use ratatui::crossterm::event::{ self, Event, KeyCode, KeyEventKind };
use ratatui::DefaultTerminal;
use ratatui::widgets::{Block, Paragraph};

mod app;
mod ui;

fn main() -> std::io::Result<()> {
    let mut terminal = ratatui::init();

    let result = run(&mut terminal);

    ratatui::restore();
    result
}


fn run(terminal: &mut DefaultTerminal) -> std::io::Result<()> {
    loop {
        terminal.draw(|frame| {
            let hello = Paragraph::new("Hello from GalaRP")
            .centered()
            .block(Block::bordered().title(" GalaRP TUI"));

            frame.render_widget(hello, frame.area());
        })?;

        if let Event::Key(key) = event::read()? {
            if key.kind == KeyEventKind::Press && key.code == KeyCode::Char('q') {
                return Ok(())
            }
        }
    }
}
