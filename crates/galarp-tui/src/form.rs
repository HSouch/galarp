//! Set up for config files.
//! Meant to preserve comments and layout
//! 
//! This needs to be done since most components are table keys, but we also have sections,
//! and [galaxy] is an array of tables.
//! 
//! Seg is one step.
//! Row is one line of the form
//!     Either a section header
//! or a field with a path.

use std::{fmt::format, path::{Path, PathBuf}};

use bigwig_io::config::Config;
use bigwig_io::snapshot::output_path;
use galarp_run::Run;
use toml_edit::{DocumentMut, Item, Table, TableLike, Value};

/// One step of a key path: a table key, or an entry of an array of tables (`[[galaxy]]`).
#[derive(Clone, Debug, PartialEq)]
pub enum Seg {
    Key(String),
    Index(usize),
}


pub fn parse_path(s: &str) -> Vec<Seg> {
    s.split('.')
        .map(|p| p.parse()
        .map_or_else(|_| Seg::Key(p.into()), Seg::Index))
        .collect()
}

pub fn path_string(path: &[Seg]) -> String {
    path.iter()
        .map(|s| match s {
            Seg::Key(k) => k.clone(),
            Seg::Index(i) => i.to_string(),
        })
        .collect::<Vec<_>>()
        .join(".")
}

#[derive(Clone, Debug)]
pub enum Row {
    Header(String),
    Field { path: Vec<Seg>, help: &'static str },
}

/// The fields every GalaRP config gets, whether or not the file sets them.
const CORE: &[(&str, &[(&str, &str)])] = &[
    (
        "run",
        &[
            ("run.output_dir", "directory for snapshots and log.csv (relative to where bigwig-tui runs)"),
            ("run.t_end_myr", "end time"),
            ("run.dt_max_myr", "log interval; a whole multiple of restricted.dt_myr"),
            ("run.snapshot_interval_myr", "snapshot interval; a whole multiple of run.dt_max_myr"),
            ("restricted.dt_myr", "RK4 step"),
            ("restricted.vertical_dispersion_km_s", "cloud sigma_z (default: the disk's sound speed)"),
            ("restricted.seed", "seed of the vertical velocities"),
        ],
    ),
    (
        "ram pressure",
        &[
            ("ram_pressure.enabled", "push the clouds with S rho |v_rel| v_rel / Sigma"),
            ("ram_pressure.surface_density_msun_pc2", "cloud surface density Sigma"),
            ("ram_pressure.density_at", "\"galaxy\" (rho at its center, as GalaRP) or \"particle\""),
            ("ram_pressure.shadow.enabled", "GalaRP's dynamic shadow"),
            ("ram_pressure.shadow.drag_coefficient", "C_D in tau = C_D A / bin_width^2"),
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

const GALAXY: &[(&str, &str)] = &[
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

pub struct ConfigDoc {
    pub path: PathBuf,
    /// The text as loaded or last saved, to mark what changed.
    pub saved: DocumentMut,
    pub doc: DocumentMut,
    pub rows: Vec<Row>,
}

impl ConfigDoc {
    pub fn load(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let doc: DocumentMut = text.parse().map_err(|e| format!("{}: {e}", path.display()))?;
        let mut c = Self { path: path.to_path_buf(), saved: doc.clone(), doc, rows: Vec::new() };
        c.build_rows();
        Ok(c)
    }

    /// A config that isn't on disk yet (`path` is where it will be saved), from `text`.
    /// Everything in it counts as an unsaved change.
    pub fn new_unsaved(path: &Path, text: &str) -> Result<Self, String> {
        let doc: DocumentMut = text.parse().map_err(|e| format!("{e}"))?;
        let mut c = Self { path: path.to_path_buf(), saved: DocumentMut::new(), doc, rows: Vec::new() };
        c.build_rows();
        Ok(c)
    }

    fn build_rows(&mut self) {
        let mut rows = Vec::new();
        for (section, fields) in CORE {
            rows.push(Row::Header(section.to_string()));
            rows.extend(fields.iter().map(|(p, help)| Row::Field { path: parse_path(p), help }));
        }
        rows.push(Row::Header("cluster gas".into()));
        let gas = get(&self.doc, &parse_path("cluster.gas")).and_then(Item::as_table_like);
        let keys: Vec<String> = gas.map(|t| t.iter().filter(|(_, v)| v.is_value()).map(|(k, _)| k.to_string()).collect()).unwrap_or_default();
        if keys.is_empty() {
            rows.push(Row::Field { path: parse_path("cluster.gas.model"), help: "\"uniform\", \"beta\", \"table\" or \"none\"" });
        }
        for k in keys {
            rows.push(Row::Field { path: vec![Seg::Key("cluster".into()), Seg::Key("gas".into()), Seg::Key(k)], help: "[cluster.gas]; see MANUAL.md for the model's keys" });
        }
        let n_gal = self.doc.get("galaxy").and_then(Item::as_array_of_tables).map_or(0, |a| a.len());
        for g in 0..n_gal {
            rows.push(
                Row::Header(
                    if g == 0 { "galaxy 0 (parent)".into() } 
                    else { "galaxy".to_string() +  &g.to_string() }
                ));
            for (p, help) in GALAXY {
                rows.push(Row::Field { path: parse_path(&format!("galaxy.{g}.{p}")), help });
            }
        }
        self.rows = rows;
    }

    pub fn name(&self) -> String {
        self.path.file_name().map_or_else(|| self.path.display().to_string(), |n| n.to_string_lossy().into())
    }

    pub fn text(&self) -> String {
        self.doc.to_string()
    }

    pub fn dirty(&self) -> bool {
        self.doc.to_string() != self.saved.to_string()
    }

    /// The value at `path` as TOML, or None if the file leaves it at its default.
    pub fn value(&self, path: &[Seg]) -> Option<String> {
        get(&self.doc, path).and_then(Item::as_value).map(show)
    }

    /// The value at `path` if it is a string, without quotes.
    pub fn str_value(&self, path: &[Seg]) -> Option<String> {
        get(&self.doc, path).and_then(Item::as_str).map(String::from)
    }

    pub fn changed(&self, path: &[Seg]) -> bool {
        get(&self.doc, path).and_then(Item::as_value).map(show) != get(&self.saved, path).and_then(Item::as_value).map(show)
    }

    /// Sets `path` from what was typed: any TOML value (`1000`, `true`, `[0, 0, 1000]`,
    /// `"beta"`), or a bare word (no quotes, brackets or spaces) taken as a string. Empty
    /// text removes the key, which puts it back to its default.
    pub fn set(&mut self, path: &[Seg], input: &str) -> Result<(), String> {
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
        if path.len() >= 2 && path[..2] == [Seg::Key("cluster".into()), Seg::Key("gas".into())] {
            self.build_rows();
        }
        Ok(())
    }

    pub fn toggle(&mut self, path: &[Seg]) -> Result<bool, String> {
        let Some(b) = get(&self.doc, path).and_then(Item::as_value).and_then(Value::as_bool).or_else(|| default_bool(path)) else {
            return Ok(false);
        };
        set(&mut self.doc, path, Some(Value::from(!b)))?;
        Ok(true)
    }

    pub fn save_as(&mut self, path: &Path) -> Result<(), String> {
        std::fs::write(path, self.doc.to_string()).map_err(|e| format!("{}: {e}", path.display()))?;
        self.path = path.to_path_buf();
        self.saved = self.doc.clone();
        Ok(())
    }

    /// Parses and validates the edited text and sets up the run (`--dry-run`).
    pub fn validate(&self) -> Validation {
        let base = self.path.parent().unwrap_or(Path::new("."));
        let text = self.text();
        let cfg = Config::from_str_in(&text, base);
        let output_dir = cfg.as_ref().map_or_else(|_| PathBuf::new(), |c| c.run.output_dir.clone());
        let result = cfg.map_err(|e| e.0).and_then(|cfg| Run::new(cfg, text).map(|r| r.summary).map_err(|e| e.to_string()));
        Validation { output_has_snapshots: !output_dir.as_os_str().is_empty() && output_path(&output_dir, 0).exists(), output_dir, result }
    }
}

/// What the edited config would run, or why it can't.
pub struct Validation {
    /// The dry-run summary, or the error.
    pub result: Result<Vec<String>, String>,
    /// `run.output_dir` (empty if the config doesn't parse), and whether a run there would
    /// overwrite snapshots.
    pub output_dir: PathBuf,
    pub output_has_snapshots: bool,
}


/// Booleans whose default is known, so a toggle works before the file sets them.
fn default_bool(path: &[Seg]) -> Option<bool> {
    matches!(path_string(path).as_str(), "ram_pressure.enabled" | "ram_pressure.shadow.enabled").then_some(false)
}

fn show(v: &Value) -> String {
    let mut v = v.clone();
    v.decor_mut().clear();
    v.to_string()
}

pub fn get<'a>(doc: &'a DocumentMut, path: &[Seg]) -> Option<&'a Item> {
    let (last, parents) = path.split_last()?;
    let mut table: &dyn TableLike = doc.as_table();
    let mut i = 0;
    while i < parents.len() {
        let Seg::Key(k) = &parents[i] else { return None };
        let item = table.get(k)?;
        table = match parents.get(i + 1) {
            Some(Seg::Index(n)) => {
                i += 1;
                item.as_array_of_tables()?.get(*n)?
            }
            _ => item.as_table_like()?,
        };
        i += 1;
    }
    match last {
        Seg::Key(k) => table.get(k),
        Seg::Index(_) => None,
    }
}

fn set(doc: &mut DocumentMut, path: &[Seg], value: Option<Value>) -> Result<(), String> {
    let bad = || format!("can't set {}", path_string(path));
    let (last, parents) = path.split_last().ok_or_else(bad)?;
    let Seg::Key(last) = last else { return Err(bad()) };
    let mut table: &mut dyn TableLike = doc.as_table_mut();
    let mut i = 0;
    while i < parents.len() {
        let Seg::Key(k) = &parents[i] else { return Err(bad()) };
        if table.get(k).is_none() {
            if value.is_none() {
                return Ok(());
            }
            table.insert(k, Item::Table(Table::new()));
        }
        let item = table.get_mut(k).ok_or_else(bad)?;
        table = match parents.get(i + 1) {
            Some(Seg::Index(n)) => {
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