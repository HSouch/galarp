

use ratatui::crossterm::event::{ self, Event, KeyCode, KeyEventKind };
use ratatui::widgets::{Block, Paragraph};
use ratatui::DefaultTerminal;

fn main() -> std::io::Result<()> {
    let mut terminal = ratatui::init();
    let res = run(&mut terminal);

    ratatui::restore();
    res
}


fn run(terminal: &mut DefaultTerminal) -> std::io::Result<()> {
    loop {
        terminal.draw(|frame| {
            let hello = Paragraph::new("Let's run some GalaRP stuff. Press q to quit.")
            .centered().block(Block::bordered().title("GalaRP"));
            frame.render_widget(hello, frame.area());
        })?;
        
        if let Event::Key(key) = event::read()? {
            if key.kind == KeyEventKind::Press && key.code == KeyCode::Char('q') {
                return Ok(())
            }
        }
    }
}