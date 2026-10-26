
use ratatui::layout::{Constraint, Flex, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span, Text};
use ratatui::widgets::{Block, Clear, List, Paragraph, Tabs, Wrap};
use ratatui::Frame;

use crate::app::{App, Mode, Tab};


pub fn draw(frame: &mut Frame, app: &mut App) {
    let [tabs, body, status] = Layout::vertical(
        [Constraint::Length(1), Constraint::Min(0), Constraint::Length(1)]).areas(frame.area());

    let titles = Tab::ALL.iter()
        .map(|t| format!("{}", t.title()));

    let selected = Tab::ALL.iter()
        .position(|&t| t == app.tab).unwrap();

    frame.render_widget(Tabs::new(titles).select(selected)
        .highlight_style(Style::new().reversed()).divider("").padding(" "," "), tabs);

    
    match app.tab {
        Tab::Configs => draw_configs(frame, body, app),
        Tab::Edit => draw_edit(frame, body, app),
        _=> placeholder(frame, body, "not built yet"),
    }

    // placeholder(frame, body, &format!("{}: not built yet", app.tab.title()));

    // frame.render_widget(Line::from("Tab switch q quit").dim(), status);
    draw_status(frame, status, app);

    match &app.mode {
        Mode::Normal => {},
        Mode::Help => draw_help(frame),
    }
}


fn placeholder(frame: &mut Frame, area: Rect, msg: &str) {
    frame.render_widget(Paragraph::new(msg.to_string())
        .dim().centered().block(Block::bordered()), area);
}   


fn draw_status(frame: &mut Frame, area: Rect, app: &App) {
    let keys = match app.tab {
        Tab::Configs => "Arrow Keys + Vim Motions move   r reload",
        _ => "",
    };
    let line = Line::from(vec![
        Span::raw(format!(" {keys}   ")), Span::raw("Tab switch   ? help   q quit").dim()
    ]);
    frame.render_widget(line, area);
}


fn draw_configs(frame: &mut Frame, area: Rect, app: &mut App) {
    let [left, right] = Layout::horizontal([
        Constraint::Percentage(40),
        Constraint::Percentage(60)
    ]).areas(area);
    
    // First do the left side with a list of available configs ///////////////////
    let items: Vec<String> = app.files.iter()
        .map(|p| app.display_path(p).into_owned()).collect();
    let title = format!(" Config Files in {}", app.cwd.display());

    let list = List::new(items)
        .block(Block::bordered().title(title))
        .highlight_style(Style::new().reversed())
        .highlight_symbol(">");

    frame.render_stateful_widget(list, left, &mut app.list);
    // ///////////////////////////////////////////////////////////////////////////


    // Now do the other side with a preview pane and config file viewer //////////
    let [summary, text] = Layout::vertical(
        [Constraint::Length(9), Constraint::Min(0)])
        .areas(right);

    let Some(check) = app.selected_check() else {
        placeholder(frame, right, "no *.toml files here or in configs/");
        return;
    };

    let (lines, color) = match &check.summary {
        Ok(lines) => (lines.iter().map(|l| Line::raw(l.clone())).collect(), Color::Green),
        Err(e) => (vec![Line::raw(e.clone())], Color::Red),
    };

    let summary_block = Block::bordered().title(" Dry Run ")
        .border_style(Style::new().fg(color));

    frame.render_widget(
        Paragraph::new(lines).block(summary_block).wrap(Wrap { trim: false}), 
        summary);

    frame.render_widget(
        Paragraph::new(check.text.clone()).block(Block::bordered().title(" File "))
            .wrap(Wrap { trim: false}), 
        text);
    // Now do the other side with a preview pane and config file viewer //////////
}   


fn draw_edit(frame: &mut Frame, area: Rect, app: &mut App) {
    
}


fn draw_help(frame: &mut Frame) {
    let text = Text::from(vec![
        Line::from("GalaRP TUI".add_modifier(Modifier::BOLD)),
        Line::raw(""),
        Line::raw("Tab / Shift-Tab   switch tabs"),
        Line::raw("↑↓ or j k         move"),
        Line::raw(""),
        Line::from("Configs".bold()),
        Line::raw("r                 reload the config list"),
        Line::raw(""),
        Line::raw("?                 this help"),
        Line::raw("q                 quit"),
        Line::raw(""),
        Line::from("any key closes this".dim()),
    ]);
    let area = popup(frame, 48, text.height() as u16 + 2);
    frame.render_widget(Paragraph::new(text).block(Block::bordered().title(" Help ")), area);

}


fn popup(frame: &mut Frame, width: u16, height: u16) -> Rect {
    let [area] = Layout::horizontal([Constraint::Length(width)]).flex(Flex::Center).areas(frame.area());
    let [area] = Layout::vertical([Constraint::Length(height)]).flex(Flex::Center).areas(area);

    frame.render_widget(Clear, area);
    area
}
