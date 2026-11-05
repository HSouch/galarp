
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};


#[derive(Default)]
pub struct LineInput {
    pub text: String,
    
    /// Cursor position in characters
    pub cursor: usize,
}


impl LineInput {
    pub fn new(text: String) -> Self {
        let cursor = text.chars().count();
        Self {
            text: text,
            cursor: cursor,
        }
    }

    fn byte(&self, i: usize) -> usize {
        self.text
            .char_indices()
            .nth(i)
            .map_or(self.text.len(), |(b, _)| b)
    }


    pub fn on_key(&mut self, key: KeyEvent) -> bool {

        let len = self.text.chars().count();

        match key.code {
            // ctrl-u erases the whole field
            KeyCode::Char('u') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                self.text.clear();
                self.cursor = 0;
            },
            // Any letter/number character writes to text
            KeyCode::Char(c) => {
                let b = self.byte(self.cursor);
                self.text.insert(b, c);
                self.cursor += 1;
            },

            KeyCode::Backspace if self.cursor > 0 => {
                self.cursor -= 1;
                let b = self.byte(self.cursor);
                self.text.remove(b);
            },

            KeyCode::Delete if self.cursor < len => {
                let b = self.byte(self.cursor);
                self.text.remove(b);
            },
            
            KeyCode::Left => self.cursor = self.cursor.saturating_sub(1),
            KeyCode::Right => self.cursor = {self.cursor + 1}.min(len),
            KeyCode::Home => self.cursor = 0,
            KeyCode::End => self.cursor = len,
            
            KeyCode::Backspace | KeyCode::Delete => {},
            _ => return false,
        }
        
        true
    }

}



