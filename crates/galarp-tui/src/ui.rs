//! UI drawing
//! 


use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, List, Paragraph, Tabs};
use ratatui::Frame;

use crate::app::{App, Tab};


pub fn draw(frame: &mut Frame, app: &mut App) {
    let [tabs, body, status] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(0),
        Constraint::Length(1)
    ]).areas(frame.area());


    let titles = Tab::ALL.iter().map(|t| format!(" {} ", t.title()));
    let selected = Tab::ALL.iter().position(|&t| t == app.tab).unwrap();

    let tab_bar = Tabs::new(titles)
        .select(selected)
        .highlight_style(Style::new().reversed())
        .divider("")
        .padding("", "");
    
    frame.render_widget(tab_bar, tabs);

    match app.tab {
        Tab::Configs => draw_configs(frame, body, app),
        _ => placeholder(frame, body, "not built yet"),
    }

    draw_status(frame, status,app);
}


fn draw_configs(frame: &mut Frame, area: Rect, app: &mut App) {
    let [left, right] = Layout::horizontal([
        Constraint::Percentage(40),
        Constraint::Percentage(60)
    ]).areas(area);

    let items: Vec<String> = app.files.iter().map(
        |p| app.display_path(p).into_owned()).collect();
    
    let title = format!(" configs in {}", app.cwd.display());

    let list = List::new(items)
        .block(Block::bordered().title(title))
        .highlight_style(Style::new().reversed())

        .highlight_symbol("> ");

    frame.render_stateful_widget(list, left, &mut app.list);

    let Some(path) = app.selected() else {
        placeholder(frame, right, "no *.toml found here or in configs/");
        return;
    };

    let text = std::fs::read_to_string(path).unwrap_or_else(|e| format!("can't read: {e}"));

    frame.render_widget(Paragraph::new(text).block(Block::bordered().title(" file ")), right);


}


fn draw_status(frame: &mut Frame, area: Rect, app: &App) {
    let keys = match app.tab {
        Tab::Configs => "Arrow-keys move   r reload",
        _ => "",
    };
    let line = Line::from(vec![Span::raw(format!(" {keys}   ")), 
        Span::raw("Tab switch   q quit").dim()]);
    
    frame.render_widget(line, area);
}


fn placeholder(frame: &mut Frame, area: Rect, msg: &str) {
    let text = Paragraph::new(msg.to_string()).dim().centered().block(Block::bordered());
    frame.render_widget(text, area);
}