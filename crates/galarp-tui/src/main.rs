

use std::time::Duration;

use ratatui::crossterm::event::{ self, Event, KeyEventKind };
use ratatui::DefaultTerminal;

mod app;            // App functionality and general keystroke handling
mod configedit;     // Editing and writing config files.
mod input;          // Line input handling
mod plotting;       // Plotting routines
mod preview;        // Getting GalaRP previews of infall (very low cost)
mod sim_runner;     // Running simulations on separate thread
mod sim_viewer;     // Viewing simulations
mod tabs;           // Individual galarp-tui tab crates
mod ui;             // UI drawing 


fn main() -> std::io::Result<()> {
    let mut app = app::App::new(std::env::current_dir()?);

    let mut terminal = ratatui::init();

    let result = run(&mut terminal, &mut app);

    ratatui::restore();
    result
}


fn run(terminal: &mut DefaultTerminal, app: &mut app::App) -> std::io::Result<()> {
    while !app.quit {
        app.tick();

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
