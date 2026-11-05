//! UI drawing functions and utilities
//! 


use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, List, Paragraph, Tabs, Wrap};
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

    
    let [summary, text] = Layout::vertical([
        Constraint::Length(9),
        Constraint::Min(0),
    ]).areas(right);

    let Some(check) = app.selected_check() else {
        placeholder(frame, right, "no .toml here or in configs/");
        return;
    };
    
    let (lines, color) = match &check.summary {
        Ok(lines) => (lines.iter().map(|l| Line::raw(l.clone())).collect(), Color::Green),
        Err(e) => (vec![Line::raw(e.clone())], Color::Red)
    };

    let summary_block = Block::bordered().title(" Run Summary ")
        .border_style(Style::new().fg(color));

    frame.render_widget(
        Paragraph::new(lines).block(summary_block).wrap(Wrap{trim: false}),  
        summary);

    frame.render_widget(
        Paragraph::new(check.text.clone()).block(Block::bordered().title(" Config File ")), 
        text);



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