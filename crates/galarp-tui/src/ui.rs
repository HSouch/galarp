//! UI drawing functions and utilities
//! 


use ratatui::layout::{Constraint, Flex, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span, Text};
use ratatui::widgets::{Block, Clear, List, Paragraph, Tabs, Wrap};
use ratatui::Frame;

use crate::app::{App, Mode, Tab};
use crate::configedit::{path_string, Row};
use crate::plotting;
use crate::plotting::timeseries::{self, YAxis};
use crate::preview::Preview;


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

    match &app.mode {
        Mode::Normal => {},
        Mode::Help => draw_help_popup(frame),

        Mode::EditField { path, input } => {
            draw_input(frame, &path_string(path), &input.text, input.cursor)
        },
        
        Mode::SaveAs { input } => {
            let title = format!("save as (relative to {})", app.cwd.display());
            draw_input(frame, &title, &input.text, input.cursor);
        },

        Mode::Confirm { question, .. } => {
            draw_confirm(frame, question);
        },

    }

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

    
    
    let [summary, text, plots] = Layout::vertical([
        Constraint::Length(10),
        Constraint::Min(0),
        Constraint::Ratio(1, 3),
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

    match &v.preview {
        Some(Ok(p)) => draw_preview(frame, plots, p, &app.projection, None),
        Some(Err(e)) => placeholder(frame, plots, &format!("no preview: {e}")),
        None => placeholder(frame, text, "no preview until the config parses"),
    }
}


fn draw_preview(frame: &mut Frame, area: Rect, p: &Preview, projection: &plotting::orbits::Projection,
    current: Option<usize>) {
    let [p_ram, phi, orbit] = Layout::horizontal([
        Constraint::Percentage(33), Constraint::Percentage(33), Constraint::Percentage(33)]
    ).areas(area);

    let t_end = p.t_myr.last().copied().unwrap_or(1.0);

    let cursor = current.and_then(|i| p.t_myr.get(i).copied());

    let series = |y: &[f64]| -> Vec<(f64, f64)> {
        p.t_myr.iter().copied().zip(y.iter().copied()).collect() 
    };

    let p_title = if p.ram_pressure_on { "P_ram [dyn/cm^2]" } else {"P_ram, off"};
    
    plotting::timeseries::time_xy_chart(frame, p_ram, p_title, &series(&p.p_ram_dyn_cm2), 
        t_end, YAxis::Auto, cursor);
    plotting::timeseries::time_xy_chart(frame, phi, "Phi_dw [deg]", &series(&p.phi_dw_deg), 
        t_end, YAxis::Auto, cursor);

    plotting::orbits::orbit_canvas(frame, orbit, &p.pos, p.r_vir_kpc, *projection, 
        current);
}


fn draw_status(frame: &mut Frame, area: Rect, app: &App) {
    let line = if let Some((msg, is_error)) = &app.message {
        Line::from(format!(" {msg}").fg(if *is_error { Color::Red } else { Color::Green }))
    } else {
        let keys = match (&app.mode, app.tab) {
            (Mode::EditField { .. }, _) => "Enter set  Esc cancel  (empty: back to the default)",
            (Mode::SaveAs { .. }, _) => "Enter save  Esc cancel",
            (Mode::Confirm { .. }, _) => "y yes  n no",
            (_, Tab::Configs) => "↑↓/jk move  Enter open  r reload",
            (_, Tab::Edit) => "↑↓/jk move  Enter edit/toggle  d default  w save",
            _ => "",
        };
        Line::from(vec![
            Span::raw(format!(" {keys}  ")), 
            Span::raw("Tab switch  ? help  q quit").dim()
        ])
    };
    frame.render_widget(line, area);
}


fn draw_input(frame: &mut Frame, title: &str, text: &str, cursor: usize) {
    let width = (frame.area().width.saturating_sub(4)).min(80);

    let area = popup(frame, width, 3);

    let block = Block::bordered().title(format!(" {title} ")).border_style(Style::new().fg(Color::Cyan));
    let inner = block.inner(area);

    let scroll = (cursor as u16).saturating_sub(inner.width.saturating_sub(1));

    frame.render_widget(Paragraph::new(text.to_string()).scroll((0, scroll)).block(block), area);
    frame.set_cursor_position((inner.x + cursor as u16 - scroll, inner.y));
}


fn draw_confirm(frame: &mut Frame, question: &str) {
    let area = popup(frame, (question.len() as u16 + 4).max(25), 4);
    let text = vec![Line::raw(question.to_string()), Line::from("y yes   n no".dim())];

    let block = Block::bordered().border_style(Style::new().fg(Color::LightBlue));
    frame.render_widget(Paragraph::new(text).block(block), area);
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


fn placeholder(frame: &mut Frame, area: Rect, msg: &str) {
    let text = Paragraph::new(msg.to_string()).dim().centered().block(Block::bordered());
    frame.render_widget(text, area);
}