//! UI drawing functions and utilities
//! 


use ratatui::layout::{Constraint, Flex, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span, Text};
use ratatui::widgets::{Block, Clear, List, Paragraph, Tabs, Wrap};
use ratatui::Frame;

use crate::app::{App, Mode, Tab};
use crate::configedit::{path_string, Row};


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
        Tab::Edit => draw_edit(frame, body, app),
       
        _ => placeholder(frame, body, "not built yet"),
    }

    draw_status(frame, status,app);

    match app.mode {
        Mode::Normal => {},
        Mode::Help => draw_help_popup(frame),
    }

}


fn popup(frame: &mut Frame, width: u16, height: u16) -> Rect {
    let [area] = Layout::horizontal(
        [Constraint::Length(width)]).flex(Flex::Center).areas(frame.area());
    let [area] = Layout::vertical([Constraint::Length(height)]).flex(Flex::Center).areas(area);
    
    // First render everything in the area clear
    clear_area(frame, area);
    area
}


fn clear_area(frame: &mut Frame, area: Rect) {
    frame.render_widget(Clear, area);
}


fn draw_help_popup(frame: &mut Frame) {
    let text = Text::from(vec![
        Line::from("galarp-tui".add_modifier(Modifier::BOLD)),
        Line::raw(""),
        Line::raw("Tab / Shift-Tab switch tabs"),
        Line::raw("↑↓ or j k move"),
        Line::raw("g / G first / last"),
        Line::raw(""),
        Line::from("Configs".bold()),
        Line::raw("r reload the config list"),
        Line::raw(""),
        Line::raw("? this help"),
        Line::raw("q quit"),
        Line::raw(""),
        Line::from("any key closes this".dim()),]);
        let area = popup(frame, 48, text.height() as u16 + 2);

    frame.render_widget(
        Paragraph::new(text).block(Block::bordered().title(" Help and Commands ")),
        area);
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


fn draw_edit(frame: &mut Frame, area: Rect, app: &mut App) {
    let Some(editor) = app.editor.as_mut() else {
        placeholder(frame, area, "Open config file from Configs tab");
        return;
    };

    let [left, right] = Layout::horizontal([
        Constraint::Percentage(50),
        Constraint::Percentage(50)
    ]).areas(area);

    let [form, help] = Layout::vertical([
            Constraint::Min(0),
            Constraint::Length(5)
    ]).areas(left);
    
    
    let width = editor.doc.rows.iter()
        .filter_map(|r| if let Row::Field { path, .. } = r 
            { Some(path_string(path).len())} else {None})
        .max()
        .unwrap_or(0);
    
    let items: Vec<Line> = editor.doc.rows.iter()
        .map(|row| match row {
            Row::Header(h) => Line::from(format!("-- {h} ")).bold(),
            Row::Field { path, .. } => {
                let key = format!("{:width$}  ", path_string(path));
                let (value, style) = match editor.doc.value(path) {
                    None => ("default".to_string(), Style::new().dim()),
                    Some(v) if editor.doc.changed(path) => (format!("{v} *"), Style::new().fg(Color::Yellow)),
                    Some(v) => (v, Style::new()),
                };
                Line::from(vec![Span::raw(key), Span::styled(value, style)])
            }
        })
        .collect();

    let title = format!(" {}{} ", editor.doc.name(), if editor.doc.not_pristine() { " [modified] "} else { ""});

    let list = List::new(items).block(Block::bordered().title(title))
        .highlight_style(Style::new().reversed());

    frame.render_stateful_widget(list, form, &mut editor.list);

    let help_text = match editor.list.selected().and_then(|i| editor.doc.rows.get(i)) {
        Some(Row::Field{ help, .. }) => help.to_string(),
        _ => String::new(),
    };
    frame.render_widget(Paragraph::new(help_text).dim()
        .wrap(Wrap { trim: true}).block(Block::bordered()), 
        help);

    
    
    let [summary, text] = Layout::vertical([
        Constraint::Length(10),
        Constraint::Min(0),
    ]).areas(right);

    let v = &editor.validation;
    let (mut lines, color): (Vec<Line>, _) = match &v.result {
        Ok(lines) => (lines.iter().map(|l| Line::raw(l.clone())).collect(), Color::Green),
        Err(e) => (vec![Line::raw(e.clone())], Color::Red)
    };

    if !v.output_dir.as_os_str().is_empty() {
        lines.push(Line::raw(""));
        lines.push(Line::raw(format!("output: {}", v.output_dir.display())));
        if v.output_has_snapshots {
            lines.push(Line::from("  already has snapshots. A run would overwrite them.".fg(Color::Yellow)));
        }
    }

    let summary_block = Block::bordered().title("  Dry Run  ").border_style(Style::new().fg(color));
    frame.render_widget(Paragraph::new(lines).wrap(Wrap {trim: false}).block(summary_block), summary);
    
    frame.render_widget(Paragraph::new(editor.doc.text()).block(Block::bordered().title(" toml")), text);
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