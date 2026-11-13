//! The data behind the TUI's config editor tab.
//!
//! The editor shows a config file as a flat list of rows: section headers
//! ("run", "ram pressure", "galaxy 0 (parent)", ...) with one editable field
//! under each. 
//! 
//! This module loads the file, decides which rows to show, and reads / writes the values.
//!
//! 1. `ConfigDoc::load` parses the file with `toml_edit`, not `toml`.
//!    `toml_edit` keeps the user's comments, spacing and key order, so saving
//!    rewrites only the values that were edited.
//! 2. `ConfigDoc` holds two copies of the document: `saved` (what is on disk)
//!    and `doc` (with the edits).
//! 3. The rows come from the tables at the bottom of this file (`CORE_FIELDS`
//!    and `GALAXY_FIELDS`), plus one row per key the file has in `[cluster.gas]`
//!    (its keys depend on the gas model, so they are read from the file).
//! 4. Each field is named by a `Vec<Segment>` path such as
//!    `galaxy.0.disk.n_particles`. `get_item` and `set_item` follow a path
//!    through the document.
//! 5. `ConfigDoc::validate` runs the edited text through the same parsing and
//!    dry run as `galarp-run`, so mistakes show while editing.
//!
//! To show a new parameter in the editor, add a `(path, help)` line to
//! `CORE_FIELDS` (for a top-level section) or `GALAXY_FIELDS` (for each
//! `[[galaxy]]`).

use std::fmt::Display;
use std::path::{Path, PathBuf};

use galarp_io::config::Config;
use galarp_io::snapshot::output_path;

use toml_edit::{DocumentMut, Item, Table, TableLike, Value};

use crate::preview::Preview;

// ---------------------------------------------------------------------------
// Paths through the document
// ---------------------------------------------------------------------------

/// One step of a path through a TOML document.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Segment {
    /// A key of a table, e.g. `run` or `t_end_myr`.
    Key(String),
    /// An entry of an array of tables, e.g. the `0` in `galaxy.0` is the first `[[galaxy]]`.
    Index(usize),
}

/// Splits a dotted path such as `"galaxy.0.disk.n_particles"` into segments.
///
/// A part that is a whole number becomes `Segment::Index` and anything else becomes
/// `Segment::Key`, so `"galaxy.0.disk"` is `[Key("galaxy"), Index(0), Key("disk")]`.
/// Keeping indices separate means `get_item` and `set_item` never have to guess
/// whether `0` is a key or an entry of `[[galaxy]]`.
pub fn parse_path(s: &str) -> Vec<Segment> {
    let mut path: Vec<Segment> = Vec::new();
    for part in s.split('.') {
        let segment: Segment = match part.parse::<usize>() {
            Ok(index) => Segment::Index(index),
            Err(_) => Segment::Key(part.to_string()),
        };
        path.push(segment);
    }
    path
}

/// The reverse of `parse_path`: `[Key("galaxy"), Index(0)]` becomes `"galaxy.0"`.
/// This is the text shown in the editor's left column.
pub fn path_string(path: &[Segment]) -> String {
    let parts: Vec<String> = path
        .iter()
        .map(|segment: &Segment| match segment {
            Segment::Key(key) => key.clone(),
            Segment::Index(index) => index.to_string(),
        })
        .collect();
    parts.join(".")
}

/// True if `path` is inside `[cluster.gas]`.
fn is_cluster_gas(path: &[Segment]) -> bool {
    let cluster_gas: [Segment; 2] = [Segment::Key("cluster".into()), Segment::Key("gas".into())];
    path.len() >= 2 && path[..2] == cluster_gas
}

// ---------------------------------------------------------------------------
// Rows of the editor list
// ---------------------------------------------------------------------------

/// One line of the editor list.
#[derive(Clone, Debug)]
pub enum Row {
    /// A section title, e.g. "run" or "galaxy 1". It can't be selected.
    Header(String),
    /// An editable value: where it is in the document, and the help text shown under the list.
    Field { path: Vec<Segment>, help: &'static str },
}

/// Builds the editor list for `doc`: the core sections, then `[cluster.gas]`,
/// then one block per `[[galaxy]]`.
fn build_rows(doc: &DocumentMut) -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();
    rows.extend(core_rows());
    rows.extend(cluster_gas_rows(doc));
    rows.extend(galaxy_rows(doc));
    rows
}

/// A header for each section of `CORE_FIELDS`, followed by its fields.
fn core_rows() -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();
    for (section, fields) in CORE_FIELDS {
        rows.push(Row::Header(section.to_string()));
        for (path, help) in fields.iter() {
            rows.push(Row::Field { path: parse_path(path), help });
        }
    }
    rows
}

/// The `[cluster.gas]` block. Which keys exist depends on the gas model
/// (`uniform`, `beta`, ...), so instead of a fixed list this shows every plain value
/// the file has in `[cluster.gas]`. Sub-tables are skipped.
///
/// If the file has no `[cluster.gas]` values at all, it shows `cluster.gas.model`
/// so there is something to edit. Setting the model adds keys, and
/// `ConfigDoc::set` rebuilds the rows when that happens.
fn cluster_gas_rows(doc: &DocumentMut) -> Vec<Row> {
    let mut rows: Vec<Row> = vec![Row::Header("cluster gas".into())];

    let gas_table: Option<&dyn TableLike> =
        get_item(doc, &parse_path("cluster.gas")).and_then(Item::as_table_like);

    let mut keys: Vec<String> = Vec::new();
    if let Some(table) = gas_table {
        for (key, item) in table.iter() {
            if item.is_value() {
                keys.push(key.to_string());
            }
        }
    }

    if keys.is_empty() {
        rows.push(Row::Field {
            path: parse_path("cluster.gas.model"),
            help: "\"uniform\", \"beta\", \"table\" or \"none\"",
        });
    }

    for key in keys {
        let path: Vec<Segment> =
            vec![Segment::Key("cluster".into()), Segment::Key("gas".into()), Segment::Key(key)];
        rows.push(Row::Field { path, help: "[cluster.gas]; see the manual for the model's keys" });
    }
    rows
}

/// One block of `GALAXY_FIELDS` per `[[galaxy]]` in the file, with paths
/// `galaxy.<index>.<field>`. Galaxy 0 is the parent.
fn galaxy_rows(doc: &DocumentMut) -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();

    let n_galaxies: usize = doc
        .get("galaxy")
        .and_then(Item::as_array_of_tables)
        .map_or(0, |galaxies| galaxies.len());

    for g in 0..n_galaxies {
        let header: String = if g == 0 {
            "galaxy 0 (parent)".to_string()
        } else {
            "galaxy ".to_string() + &g.to_string()
        };
        rows.push(Row::Header(header));

        for (field, help) in GALAXY_FIELDS {
            // e.g. galaxy.1.disk.mass_msun
            let mut path: Vec<Segment> = vec![Segment::Key("galaxy".into()), Segment::Index(g)];
            path.extend(parse_path(field));
            rows.push(Row::Field { path, help });
        }
    }
    rows
}

// ---------------------------------------------------------------------------
// The document being edited
// ---------------------------------------------------------------------------

/// A config file open in the editor.
pub struct ConfigDoc {
    /// Where the file is. Changes on "save as".
    pub path: PathBuf,
    /// The document as it is on disk (as loaded, or as last saved).
    pub saved: DocumentMut,
    /// The document with the user's edits.
    pub doc: DocumentMut,
    /// The editor list, built by `build_rows`.
    pub rows: Vec<Row>,
}

impl ConfigDoc {
    /// Reads and parses the file at `path`.
    pub fn load(path: &Path) -> Result<Self, String> {
        let text: String = std::fs::read_to_string(path).map_err(|e| with_path(path, e))?;
        let doc: DocumentMut = text.parse().map_err(|e| with_path(path, e))?;

        let rows: Vec<Row> = build_rows(&doc);
        Ok(Self { path: path.to_path_buf(), saved: doc.clone(), doc, rows })
    }

    /// The file name for the editor's title, e.g. `face_on.toml`.
    pub fn name(&self) -> String {
        match self.path.file_name() {
            Some(name) => name.to_string_lossy().into(),
            None => self.path.display().to_string(),
        }
    }

    /// The edited document as TOML text.
    pub fn text(&self) -> String {
        self.doc.to_string()
    }

    /// True if there are unsaved edits.
    pub fn not_pristine(&self) -> bool {
        self.doc.to_string() != self.saved.to_string()
    }

    // ----- reading values -----

    /// The value at `path` as TOML text (strings keep their quotes),
    /// or `None` if the file doesn't set it, so the default is used.
    pub fn value(&self, path: &[Segment]) -> Option<String> {
        value_text_at(&self.doc, path)
    }

    /// The value at `path` if it is a string, without quotes.
    pub fn str_value(&self, path: &[Segment]) -> Option<String> {
        get_item(&self.doc, path).and_then(Item::as_str).map(String::from)
    }

    /// True if the value at `path` was edited since loading or saving.
    /// Adding or removing the key counts as a change.
    pub fn changed(&self, path: &[Segment]) -> bool {
        let now: Option<String> = value_text_at(&self.doc, path);
        let before: Option<String> = value_text_at(&self.saved, path);
        now != before
    }

    // ----- editing values -----

    /// Sets `path` to what the user typed, see `parse_input`.
    /// Empty input removes the key, so the default is used.
    pub fn set(&mut self, path: &[Segment], input: &str) -> Result<(), String> {
        let value: Option<Value> = parse_input(input)?;
        set_item(&mut self.doc, path, value)?;

        // The [cluster.gas] rows are read from the file (see `cluster_gas_rows`),
        // so changing a gas value, e.g. the model, can change which rows to show.
        if is_cluster_gas(path) {
            self.rows = build_rows(&self.doc);
        }
        Ok(())
    }

    /// Sets `path` to the string `s`, keeping the old value's comments.
    /// Unlike `set`, `s` is never parsed as TOML.
    pub fn set_str(&mut self, path: &[Segment], s: &str) -> Result<(), String> {
        set_item(&mut self.doc, path, Some(Value::from(s)))
    }

    /// Flips the boolean at `path`.
    ///
    /// Returns `Ok(false)` and changes nothing if the value isn't a boolean (or isn't
    /// set). The caller then opens the text input instead.
    pub fn toggle(&mut self, path: &[Segment]) -> Result<bool, String> {
        let current: Option<bool> =
            get_item(&self.doc, path).and_then(Item::as_value).and_then(Value::as_bool);

        let Some(b) = current else {
            return Ok(false);
        };

        set_item(&mut self.doc, path, Some(Value::from(!b)))?;
        Ok(true)
    }

    // ----- saving and checking -----

    /// Writes the edited document to `path`. Afterwards `path` is the doc's path
    /// and the written document is the `saved` copy.
    pub fn save_as(&mut self, path: &Path) -> Result<(), String> {
        std::fs::write(path, self.doc.to_string()).map_err(|e| with_path(path, e))?;

        self.path = path.to_path_buf();
        self.saved = self.doc.clone();
        Ok(())
    }

    /// Checks the edited text: parses it as a `Config` (relative paths are taken from
    /// the file's folder), does a dry run, and makes the preview.
    pub fn validate(&self) -> Validation {
        let base_dir: &Path = self.path.parent().unwrap_or(Path::new("."));
        let text: String = self.text();

        let config = Config::from_str_in(&text, base_dir);

        // These two need a config that parsed. If it didn't, the output dir is
        // empty and there is no preview; the error is in `result`.
        let output_dir: PathBuf = match &config {
            Ok(c) => c.run.output_dir.clone(),
            Err(_) => PathBuf::new(),
        };
        let preview: Option<Result<Preview, String>> = config.as_ref().ok().map(Preview::new);

        // The config error, or the dry run's summary (or error).
        let result: Result<Vec<String>, String> = match config {
            Err(e) => Err(e.0),
            Ok(c) => galarp_run::dry_run(&c).map_err(|e| e.to_string()),
        };

        let output_has_previous: bool =
            !output_dir.as_os_str().is_empty() && has_previous_snapshots(&output_dir);

        Validation { result, output_dir, output_has_previous, preview }
    }
}

/// The result of `ConfigDoc::validate`, shown in the editor.
pub struct Validation {
    /// The dry run's summary lines, or why the config is invalid.
    pub result: Result<Vec<String>, String>,

    /// `run.output_dir`, empty if the config didn't parse.
    pub output_dir: PathBuf,
    /// True if `output_dir` already has output, so running would overwrite it.
    pub output_has_previous: bool,

    /// The preview of the initial conditions, or `None` if the config didn't parse.
    pub preview: Option<Result<Preview, String>>,
}

/// True if `dir` has the output of an earlier run: the first snapshot,
/// `log.csv`, or the copy of the config.
pub fn has_previous_snapshots(dir: &Path) -> bool {
    output_path(dir, 0).exists() || dir.join("log.csv").exists() || dir.join("config.toml").exists()
}

// ---------------------------------------------------------------------------
// Helpers for values and errors
// ---------------------------------------------------------------------------

/// Turns what the user typed into a value. Empty input is `None` (remove the key).
///
/// The input is parsed as a TOML value first (`1.5`, `true`, `"beta"`, `[1, 2, 3]`).
/// If that fails and the input is a bare word such as `beta`, it is taken as the
/// string `"beta"` so the quotes can be left out. Anything else that isn't TOML
/// (e.g. `[1, 2`) is an error.
fn parse_input(input: &str) -> Result<Option<Value>, String> {
    let input: &str = input.trim();
    if input.is_empty() {
        return Ok(None);
    }

    match input.parse::<Value>() {
        Ok(value) => Ok(Some(value)),
        Err(_) if is_bare_word(input) => Ok(Some(Value::from(input))),
        Err(_) => Err(input.to_string() + " is not a TOML value"),
    }
}

/// True if `s` has none of the characters that make up TOML strings, arrays and
/// tables (quotes, brackets, braces, `=`, `,` or a space).
fn is_bare_word(s: &str) -> bool {
    const TOML_PUNCTUATION: &str = "\"'[]{}=, ";
    !s.chars().any(|c: char| TOML_PUNCTUATION.contains(c))
}

/// The value at `path` in `doc` as TOML text, or `None` if it's unset or a table.
fn value_text_at(doc: &DocumentMut, path: &[Segment]) -> Option<String> {
    get_item(doc, path).and_then(Item::as_value).map(value_text)
}

/// A value as TOML text, without its decor (the spaces and comments
/// `toml_edit` keeps around it), e.g. `1.5` rather than ` 1.5 # kpc`.
fn value_text(value: &Value) -> String {
    let mut bare: Value = value.clone();
    bare.decor_mut().clear();
    bare.to_string()
}

/// An error message that starts with the file, e.g. `configs/a.toml: not found`.
fn with_path(path: &Path, error: impl Display) -> String {
    path.display().to_string() + ": " + &error.to_string()
}

// ---------------------------------------------------------------------------
// Following a path through the document
// ---------------------------------------------------------------------------
//   - `Key(k)` steps into table `k`.
//   - `Key(k)` followed by `Index(n)` steps into the n-th entry of the array of
//     tables `k`, and uses up both segments.
// The last segment must be a `Key`: the value to get or set.

/// The item at `path`, or `None` if any step is missing or isn't a table.
pub fn get_item<'a>(doc: &'a DocumentMut, path: &[Segment]) -> Option<&'a Item> {
    let (last, parents) = path.split_last()?;

    let mut table: &dyn TableLike = doc.as_table();
    let mut i: usize = 0;

    while i < parents.len() {
        // A table step must start with a key (an Index here would have no array).
        let Segment::Key(key) = &parents[i] else { return None };
        let item: &Item = table.get(key)?;

        table = match parents.get(i + 1) {
            Some(Segment::Index(n)) => {
                i += 1; // the index is used up here as well
                item.as_array_of_tables()?.get(*n)?
            }
            _ => item.as_table_like()?,
        };

        i += 1;
    }

    match last {
        Segment::Key(key) => table.get(key),
        Segment::Index(_) => None,
    }
}

/// Sets (`Some`) or removes (`None`) the value at `path`.
///
/// - Missing tables on the way are created, as `[a.b]` tables. Missing arrays of
///   tables are not: setting `galaxy.3.name` with two galaxies is an error.
/// - Removing a value whose table doesn't exist does nothing.
/// - Replacing a value keeps the old value's decor, so a trailing `# comment` survives.
fn set_item(doc: &mut DocumentMut, path: &[Segment], value: Option<Value>) -> Result<(), String> {
    let cant_set = || "can't set ".to_string() + &path_string(path);

    let (last, parents) = path.split_last().ok_or_else(cant_set)?;
    let Segment::Key(last_key) = last else { return Err(cant_set()) };

    let mut table: &mut dyn TableLike = doc.as_table_mut();
    let mut i: usize = 0;

    while i < parents.len() {
        let Segment::Key(key) = &parents[i] else { return Err(cant_set()) };

        if table.get(key).is_none() {
            if value.is_none() {
                return Ok(()); // nothing to remove
            }
            table.insert(key, Item::Table(Table::new()));
        }
        let item: &mut Item = table.get_mut(key).ok_or_else(cant_set)?;

        table = match parents.get(i + 1) {
            Some(Segment::Index(n)) => {
                i += 1; // the index is used up here as well
                item.as_array_of_tables_mut().and_then(|array| array.get_mut(*n)).ok_or_else(cant_set)?
            }
            _ => item.as_table_like_mut().ok_or_else(cant_set)?,
        };
        i += 1;
    }

    match value {
        None => {
            table.remove(last_key);
        }
        Some(mut new_value) => {
            if let Some(old_value) = table.get(last_key).and_then(Item::as_value) {
                *new_value.decor_mut() = old_value.decor().clone();
            }
            table.insert(last_key, Item::Value(new_value));
        }
    }
    Ok(())
}


/// A field: its dotted path and the help text shown when it is selected.
pub type FieldSpec = (&'static str, &'static str);

/// A section: its header and its fields.
pub type SectionSpec = (&'static str, &'static [FieldSpec]);

pub const CORE_FIELDS: &[SectionSpec] = &[
    (
        "run",
        &[
            ("run.output_dir", "directory for snapshots and log.csv (relative to the configuration file location)"),
            ("run.t_end_myr", "end time"),
            ("run.dt_max_myr", "log interval; a whole multiple of galarp.dt_myr"),
            ("run.snapshot_interval_myr", "snapshot interval; a whole multiple of run.dt_max_myr"),
            ("box.side_kpc", "side of the box written to snapshot headers, centered on the parent"),
            ("galarp.dt_myr", "RK4 step"),
            ("galarp.vertical_dispersion_km_s", "cloud sigma_z (default: the disk's sound speed)"),
            ("galarp.seed", "seed of the vertical velocities"),
        ],
    ),
    (
        "ram pressure",
        &[
            ("ram_pressure.enabled", "push the clouds with S rho |v_rel| v_rel / Sigma"),
            ("ram_pressure.surface_density_msun_pc2", "cloud surface density Sigma"),
            ("ram_pressure.density_at", "\"galaxy\" (rho at its center, as GalaRP) or \"particle\""),
            ("ram_pressure.shadow.enabled", "GalaRP's dynamic shadow"),
            ("ram_pressure.drag_coefficient", "drag coefficient on clouds"),
            ("ram_pressure.shadow.bin_width_kpc", "side of the columns along the wind"),
        ],
    ),
    (
        "diagnostics",
        &[
            ("diagnostics.disk_region_radius_kpc", "disk region radius for in_disk_fraction"),
            ("diagnostics.disk_region_half_height_kpc", "disk region half height for in_disk_fraction"),
        ],
    ),
];

pub const GALAXY_FIELDS: &[FieldSpec] = &[
    ("name", "label in the log and summaries"),
    ("inclination_deg", "disk tilt from the cluster z axis (or give disk_normal)"),
    ("azimuth_deg", "azimuth of the tilt"),
    ("disk.n_particles", "number of clouds"),
    ("disk.mass_msun", "disk gas mass"),
    ("disk.scale_length_kpc", "exponential scale length"),
    ("disk.r_max_kpc", "disk truncation radius"),
    ("disk.sound_speed_km_s", "sets the default sigma_z"),
    ("orbit.mode", "\"constant_velocity\", \"integrate\" or \"table\""),
    ("orbit.position_kpc", "start position in the cluster frame"),
    ("orbit.velocity_km_s", "start velocity in the cluster frame; the wind is minus this"),
];
