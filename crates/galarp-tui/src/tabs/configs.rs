use std::path::{Path, PathBuf};


pub fn find(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();

    for d in [dir.to_path_buf(), dir.join("configs")] {

        let Ok(entries) = std::fs::read_dir(&d) else { continue };

        let mut files: Vec<PathBuf> = entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "toml") && is_config(p))
            .collect();

        files.sort();
        out.extend(files);
    }

    out.dedup();
    out
}

/// Determine whether the given file is a proper config file.
/// Every GalaRP config file needs at least one [[galaxy]] config region.
fn is_config(path: &Path) -> bool {
    std::fs::read_to_string(path)
        .is_ok_and(|t| t.lines().any(|l| l.trim_start().starts_with("[[galaxy]]")))
}
