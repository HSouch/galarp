
use std::path::{Path, PathBuf};

use galarp_io::config::Config;
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

        rows.push(Row::Header("cluster gas".into()));

        let gas = get(&self.doc, &parse_path("cluster.gas")).and_then(Item::as_table_like);
        let keys: Vec<String> = gas
            .map(|t| t.iter().filter(|(_, v)| v.is_value()).map(|(k, _)| k.to_string()).collect())
            .unwrap_or_default();
        if keys.is_empty() {
            rows.push(
                Row::Field { path: parse_path("cluster.gas.model"), 
                help: "\"uniform\", \"beta\", \"table\" or \"none\"" });
        }
        for k in keys {
            let path = vec![Segment::Key("cluster".into()), Segment::Key("gas".into()), Segment::Key(k)];
            rows.push(Row::Field { path, help: "[cluster.gas]; see the manual for the model's keys" });
        }

        // One block of GALAXY rows per [[galaxy]].
        let n_gal = self.doc.get("galaxy").and_then(
            Item::as_array_of_tables).map_or(0, |a| a.len());
        for g in 0..n_gal {
            rows.push(
                Row::Header(if g == 0 { "galaxy 0 (parent)".into() } 
                else 
                { format!("galaxy {g}").to_string() }));
            
            for (p, help) in GALAXY {
                rows.push(Row::Field { path: parse_path(&format!("galaxy.{g}.{p}")), help });
            }
        }
       
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


pub const CORE: &[(&str, &[(&str, &str)])] = &[
    (
        "run",
        &[
            ("run.output_dir", "directory for snapshots and log.csv (relative to where galarp-tui runs)"),
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

pub const GALAXY: &[(&str, &str)] = &[
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
