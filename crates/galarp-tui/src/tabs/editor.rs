
use ratatui::widgets::ListState;

use crate::form::{ConfigDoc, Row, Seg, Validation};

pub struct Editor {
    pub doc: ConfigDoc,
    pub list: ListState,
    pub validation: Validation,
}


impl Editor {
    pub fn new(doc: ConfigDoc) -> Self {
        let validation = doc.validate();
        let first = doc.rows.iter().position(|r| matches!(r, Row::Field { .. }));
        Self { doc, list: ListState::default().with_selected(first), validation}
    }

    pub fn selected_path(&self) -> Option<&[Seg]> {
        match self.doc.rows.get(self.list.selected()?)? {
            Row::Field { path, .. } => Some(path),
            Row::Header(_) => None
        }
    }

    pub fn move_by(&mut self, step: isize) {
        let rows = &self.doc.rows;
        let mut i = self.list.selected().unwrap_or(0) as isize;

        loop {
            i += step;

            if i < 0 || i > rows.len() as isize {
                return;
            }

            if matches!(rows[i as usize], Row::Field { .. }) {
                self.list.select(Some(i as usize));
                return;
            }
        }
    }

    pub fn revalidate(&mut self) {
        self.validation = self.doc.validate();
    }

}