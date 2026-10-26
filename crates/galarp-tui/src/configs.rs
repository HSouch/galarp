use std::path::{Path, PathBuf};

use galarp_io::config::Config;
use galarp_run::Run;


pub fn find(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();

    for d in [dir.to_path_buf(), dir.join("configs")] {
        let Ok(entries) = std::fs::read_dir(&d) else { continue };

        let mut files: Vec<PathBuf> = entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "toml") && is_config(p)).collect();
        files.sort();
        out.extend(files);
    }

    out.dedup();
    out
}


fn is_config(path: &Path) -> bool {
    std::fs::read_to_string(path)
        .is_ok_and( |t| t.lines()
            .any(|l| l.trim_start()
                .starts_with("[[galaxy]]"))
        )
}

pub struct Check {
    pub summary: Result<Vec<String>, String>,
    pub text: String,
}


pub fn check(path: &Path) -> Check {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| format!("cant read: {e}"));

    let base = path.parent().unwrap_or(Path::new("."));
    
    let summary = Config::from_str_in(&text, base)
        .map_err(|e| e.0).and_then(|cfg| Run::new(cfg, text.clone())
        .map(|r| r.summary).map_err(|e| e.to_string()));
    
    Check{summary, text}
}