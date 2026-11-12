
use ratatui::widgets::ListState;

use crate::configedit::{ConfigDoc, Row, Segment, Validation};


pub struct Editor {
    pub doc: ConfigDoc,

    pub list: ListState,

    /// 
    pub validation: Validation,
}


impl Editor {
    pub fn new(doc: ConfigDoc) -> Self {
        let validation = doc.validate();

        let first = doc.rows.iter().position(|r| matches!(r, Row::Field { .. }));

        Self {
            doc,
            list: ListState::default().with_selected(first),
            validation
        }
    }

    pub fn selected_path(&self) -> Option<&[Segment]> {
        match self.doc.rows.get(self.list.selected()?)? {
            Row::Field { path, .. } => Some(path),
            Row::Header(_) => None,
        }
    }

    /// Move through the editor by `step` steps.
    pub fn move_by(&mut self, step: isize) {
        let fields: Vec<usize> = self.doc.rows.iter().enumerate()
            .filter(|(_, r)| matches!(r, Row::Field { .. }))
            .map(|(i, _)| i)
            .collect();
        if fields.is_empty() { return; }

        let current = self.list.selected().unwrap_or(fields[0]);
        let k = fields.iter().position(|&i| i >= current).unwrap_or(fields.len() - 1) as isize;
        let k = (k + step).clamp(0, fields.len() as isize - 1) as usize;

        self.list.select(Some(fields[k]));
    }

    /// Reselect a path 
    /// Used to ensure selected row isn't on a header after edit.
    pub fn reselect(&mut self, path: &[Segment]) {
        let rows = &self.doc.rows;

        let is_field = |i: usize| matches!(rows[i], Row::Field { .. });
        
        let i = rows.iter()
            .position(|r| matches!(r, Row::Field{ path: p, ..} if p.as_slice() == path))
            .or_else(|| {
                let old = self.list.selected().unwrap_or(0).min(rows.len().saturating_sub(1));
                (old..rows.len()).chain((0..old).rev()).find(|&i| is_field(i))
            });

            self.list.select(i);
    }

    pub fn revalidate(&mut self) {
        self.validation = self.doc.validate();
    }
}