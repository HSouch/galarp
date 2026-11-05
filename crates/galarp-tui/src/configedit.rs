
use std::path::{Path, PathBuf};

use galarp_io::config::{Config, CORE, GALAXY};
use galarp_io::snapshot::output_path;
use galarp_run::Run;

use toml_edit::{DocumentMut, Item, Table, TableLike, Value};


/// Either a table key, or an entry of an array of tables
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Segment {
    Key(String),
    Index(usize),
}

/// Keys like galaxy.0.disk.n_particles is a path through the document.
/// Since galaxy is an array of tables, the index needs to reference an entry.
/// This way, we don't need to guess whether 0 in this case is an index or not.
pub fn parse_path(s: &str) -> Vec<Segment> {
    s.split('.').map(
        |p| p.parse().map_or_else(|_| Segment::Key(p.into()), Segment::Index)
    ).collect()
}


pub fn path_string(path: &[Segment]) -> String {
    path.iter()
        .map(|s| match s {
            Segment::Key(k) => k.clone(),
            Segment::Index(i) => i.to_string(),
        })
        .collect::<Vec<_>>()
        .join(".")
}


#[derive(Clone, Debug)]
pub enum Row {
    Header(String),
    Field { path: Vec<Segment>, help: &'static str}
}


pub struct ConfigDoc {
    pub path: PathBuf,
    /// The document that was loaded
    pub saved: DocumentMut,
    /// The document being edited
    pub doc: DocumentMut,
    pub rows: Vec<Row>,
}


impl ConfigDoc {
    pub fn load(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let doc: DocumentMut = text.parse().map_err(|e| format!("{}: {e}", path.display()))?;
        let mut c = Self { 
            path: path.to_path_buf(), 
            saved: doc.clone(), 
            doc, 
            rows:Vec::new() };
        c.build_rows();
        Ok(c)
    }

    fn build_rows(&mut self) {
        let mut rows = Vec::new();

        // For each section, add section to rows
        // Then get all of the fields for the section and extend.
        for (section, fields) in CORE {
            rows.push(Row::Header(section.to_string()));
            rows.extend(fields.iter().map(|(p, help)| Row::Field { path: parse_path(p), help }));
        }

        // This could use more work
        
        self.rows = rows;

    }

    pub fn name(&self) -> String {
        self.path.file_name().map_or_else(
            || self.path.display().to_string(), 
            |n| n.to_string_lossy().into())
    }

    pub fn text(&self) -> String {
        self.doc.to_string()
    }

    pub fn not_pristine(&self) -> bool {
        self.doc.to_string() != self.saved.to_string()
    }

    /// The value at `path` as a TOML entry, or None if the file leaves it as the default.
    pub fn value(&self, path: &[Segment]) -> Option<String> {
        get(&self.doc, path).and_then(Item::as_value).map(show)
    }

    /// Check whether value has changed between current edited version and original
    /// 
    pub fn changed(&self, path: &[Segment]) -> bool {
        let now = get(&self.doc, path).and_then(Item::as_value).map(show);
        
        let before = get(&self.saved, path).and_then(Item::as_value).map(show);

        now != before

    }

    pub fn set(&mut self, path: &[Segment], input: &str) -> Result<(), String> {
        let input = input.trim();
        let value = if input.is_empty() {
            None
        } else {
            match input.parse::<Value>() {
                Ok(v) => Some(v),
                Err(_) if !input.contains(|c: char| "\"'[]{}=, ".contains(c)) => Some(Value::from(input)),
                Err(_) => return Err(format!("{input} is not a TOML value")),
            }
        };
        set(&mut self.doc, path, value)?;
        if path.len() >= 2 && path[..2] == [Segment::Key("cluster".into()), Segment::Key("gas".into())] {
            self.build_rows();
        }
        Ok(())
    }

    /// Flip a boolean value
    pub fn toggle(&mut self, path: &[Segment]) -> Result<bool, String> {
        let current = get(&self.doc, path)
            .and_then(Item::as_value).and_then(Value::as_bool);

        let Some(b) = current else {
            return Ok(false)
        };

        set(&mut self.doc, path, Some(Value::from(!b)))?;

        Ok(true)
    }

    /// Writes document, and makes written path the doc's path.
    /// Also makes the written text the `saved` copy.
    pub fn save_as(&mut self, path: &Path) -> Result<(), String> {
        std::fs::write(path, self.doc.to_string())
            .map_err(|e| format!("{}: {e}", path.display()))?;

        self.path = path.to_path_buf();
        self.saved = self.doc.clone();
        Ok(())
    }

    pub fn validate(&self) -> Validation {
        let base = self.path.parent().unwrap_or(Path::new("."));
        let text = self.text();

        let cfg = Config::from_str_in(&text, base);

        let output_dir = cfg.as_ref().map_or_else(
            |_| PathBuf::new(),
            |c| c.run.output_dir.clone()
        );

        let result = cfg
            .map_err(|e| e.0)
            .and_then(|cfg| Run::new(cfg, text).map(|r| r.summary)
            .map_err(|e| e.to_string())
        );
        
        // Check if output already has snapshots
        let output_has_snapshots = !output_dir.as_os_str().is_empty() && 
            output_path(&output_dir, 0).exists();

        Validation { result, output_dir, output_has_snapshots }
        
    }

    /// Value at `path` if it is a string, without quotes.
    pub fn str_value(&self, path: &[Segment]) -> Option<String> {
        get(&self.doc, path).and_then(Item::as_str).map(String::from)
    }
}


fn show(v: &Value) -> String {
    let mut v = v.clone();
    v.decor_mut().clear();

    v.to_string()
}


/// Walk through the tables
/// Step through array of tables when next segment is an index
pub fn get<'a>(doc: &'a DocumentMut, path: &[Segment]) -> Option<&'a Item> {
    let (last, parents) = path.split_last()?;

    let mut table: &dyn TableLike = doc.as_table();

    let mut i = 0;

    while i < parents.len() {

        let Segment::Key(k) = &parents[i] else { return None };
        let item = table.get(k)?;

        table = match parents.get(i + 1) {
            Some(Segment::Index(n)) => {
                i += 1;
                item.as_array_of_tables()?.get(*n)?
            }
            _ => item.as_table_like()?,
        };

        i += 1;
    }

    match last {
        Segment::Key(k) => table.get(k),
        Segment::Index(_) => None,
    }
}


fn set(doc: &mut DocumentMut, path: &[Segment], value: Option<Value>) -> Result<(), String> {
    let bad = || format!("can't set {}", path_string(path));
    
    let (last, parents) = path.split_last().ok_or_else(bad)?;
    
    let Segment::Key(last) = last else { return Err(bad()) };
    
    let mut table: &mut dyn TableLike = doc.as_table_mut();
    let mut i = 0;
    
    while i < parents.len() {
        let Segment::Key(k) = &parents[i] else { return Err(bad()) };
        if table.get(k).is_none() {
            if value.is_none() {
                return Ok(());
            }
            table.insert(k, Item::Table(Table::new()));
        }
        let item = table.get_mut(k).ok_or_else(bad)?;
        table = match parents.get(i + 1) {
            Some(Segment::Index(n)) => {
                i += 1;
                item.as_array_of_tables_mut().and_then(|a| a.get_mut(*n)).ok_or_else(bad)?
            }
            _ => item.as_table_like_mut().ok_or_else(bad)?,
        };
        i += 1;
    }
    match value {
        None => {
            table.remove(last);
        }
        Some(mut v) => {
            if let Some(old) = table.get(last).and_then(Item::as_value) {
                *v.decor_mut() = old.decor().clone();
            }
            table.insert(last, Item::Value(v));
        }
    }
    Ok(())
}


pub struct Validation {
    pub result: Result<Vec<String>, String>,

    pub output_dir: PathBuf,
    pub output_has_snapshots: bool,
}


