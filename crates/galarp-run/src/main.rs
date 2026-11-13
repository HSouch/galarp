//! Run GalaRP stripping sim from a TOML config.
//! 

use std::ops::ControlFlow;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};
use std::collections::HashMap;

use std::io::{self, BufRead, Write};

use galarp_io::config::Config;
use galarp_io::snapshot::output_path;
use galarp_run::{Error, Run};
use clap::Parser;


#[derive(Parser)]
#[command(about = "Run Galarp from a TOML config")]
struct Args {
    config: PathBuf,

    #[arg(long)]
    dry_run: bool,

    #[arg(long)]
    overwrite: bool,

    #[arg(long)]
    all: bool,

    #[arg(long, short)]
    yes: bool,
}


fn main() -> ExitCode {
    let args = Args::parse();

    let is_dir = args.config.is_dir();
    if args.all && !is_dir {
        eprintln!("error: --all expects a directory, got {}", args.config.display());
        return ExitCode::FAILURE;
    } 
    if !args.all && is_dir {
        eprintln!("error: {} is a directory. Pass --all to run all config files", args.config.display());
        return ExitCode::FAILURE;
    }

    if args.all { return run_batch(&args); }

    match run_one(&args.config, args.dry_run, args.overwrite) {
        Ok(_) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Run a single instance of GalaRP
fn run_one(config: &Path, dry_run: bool, overwrite: bool) -> Result<Duration, Error> {
    let cfg = Config::from_file(config)?;
    let config_text = std::fs::read_to_string(config)?;

    let out = cfg.run.output_dir.clone();

    if !dry_run && output_path(&out, 0).exists() && !overwrite {
        return Err(
            format!("{} already has snapshots. Pass --overwrite or change run.output_dir", 
            out.display()).into() );
    }
    let t0 = Instant::now();
    println!("config: {}", config.display());

    let mut run = Run::new(cfg, config_text)?;

    for line in &run.summary {
        println!("{line}");
    }

    if dry_run {
        return Ok(t0.elapsed());
    }

    run.execute(|p| {
        if p.chunk > 0 && p.snapshot.is_some() {
            println!(
                "t = {:6.1} Myr: parent cloud mass in disk region {:.3}, {:.1} s", 
                p.t_myr, p.row[1], t0.elapsed().as_secs_f64()
            );
        }
        ControlFlow::Continue(())
    })?;

    Ok(t0.elapsed())

}


fn run_batch(args: &Args) -> ExitCode {
    let configs = match collect_configs(&args.config) {
        Ok(c) => c,
        Err(e) => { eprintln!("error: {e}"); return ExitCode::FAILURE; }
    };
    let plans = plan(&configs);

    let t_batch = Instant::now();

    let existing: Vec<&Plan> = plans.iter().filter(|p| p.has_results).collect();
    if !existing.is_empty() {
    if args.overwrite {
        println!("{} of {} configs already have results:", existing.len(), plans.len());
        for p in &existing {
            println!("  {}", p.path.display());
        }
        if !args.dry_run && !args.yes && !confirm("Overwrite them?") {
            eprintln!("aborted; nothing was run");
            return ExitCode::FAILURE;
        }
    } else {
        println!(
            "{} configs already have results and will be skipped (pass --overwrite to rerun them)", 
            existing.len());
        }
    }

    let n = plans.len();
    let mut results: Vec<(&Plan, Outcome)> = Vec::with_capacity(n);


    for (i, plan) in plans.iter().enumerate() {
        println!("\n[{}/{n}] {}", i + 1, plan.path.display());

        let outcome = match &plan.output_dir {
            Err(e) => Outcome::Failed(e.clone()),
            Ok(_) if plan.has_results && !args.overwrite => Outcome::Skipped,
            Ok(_) => match std::panic::catch_unwind(|| run_one(&plan.path, args.dry_run, args.overwrite)) {
                Ok(Ok(elapsed)) => Outcome::Ran(elapsed),
                Ok(Err(e)) => Outcome::Failed(e.to_string()),
                Err(_) => Outcome::Failed("panicked (message above)".into()),
            },
        };

        match &outcome {
            Outcome::Failed(e) => eprintln!("  failed: {e}"),
            Outcome::Skipped => println!("  skipped: already has results"),
            Outcome::Ran(_) => {}
        }
        results.push((plan, outcome));
    }

    print_summary(&results, args.dry_run, t_batch.elapsed());

    if results.iter().any(|(_, o)| matches!(o, Outcome::Failed(_))) {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}


/// Gather all toml config files inside `dir`, sorted by name.
fn collect_configs(dir: &Path) -> Result<Vec<PathBuf>, Error> {
    let mut configs: Vec<PathBuf> = std::fs::read_dir(dir)?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| p.is_file() && p.extension().is_some_and(|ext| ext == "toml"))
        .collect();
    configs.sort();
    if configs.is_empty() {
        return Err(format!("no .toml configs in {}", dir.display()).into());
    }
    Ok(configs)
}


/// What the pre-scan found for one config.
struct Plan {
    path: PathBuf,
    /// Either holds the output dir, or why this config cannot be run.
    output_dir: Result<PathBuf, String>,
    /// The output directory already contains results.
    has_results: bool,
}


fn plan(configs: &[PathBuf]) -> Vec<Plan> {
    let mut plans: Vec<Plan> = configs
        .iter()
        .map(|path| match Config::from_file(path) {
            Ok(cfg) => {
                let out = cfg.run.output_dir;
                let has_results = output_path(&out, 0).exists();
                Plan { path: path.clone(), output_dir: Ok(out), has_results }
            }
            Err(e) => Plan { path: path.clone(), output_dir: Err(e.to_string()), has_results: false },
        })
        .collect();

    // Configs sharing an output directory would overwrite each other: refuse to run any of them.
    let mut users: HashMap<PathBuf, Vec<usize>> = HashMap::new();
    for (i, p) in plans.iter().enumerate() {
        if let Ok(out) = &p.output_dir {
            users.entry(out.clone()).or_default().push(i);
        }
    }
    for (out, idx) in users.into_iter().filter(|(_, idx)| idx.len() > 1) {
        let names: Vec<String> = idx.iter()
        .map(|&i| plans[i].path.display().to_string()).collect();
        let msg = format!("shares output_dir {} with {}", out.display(), names.join(", "));
        for &i in &idx {
            plans[i].output_dir = Err(msg.clone());
        }
    }
    plans
}


/// Enum to handle outcomes of GalaRP events
/// 
enum Outcome {
    Ran(Duration),
    Skipped,
    Failed(String),
}

/// User-confirmed statement.
fn confirm(question: &str) -> bool {
    print!("{question} [y/N] ");
    let _ = io::stdout().flush();
    let mut line = String::new();
    
    if io::stdin().lock().read_line(&mut line).is_err() {
        return false;
    }

    matches!(line.trim().to_ascii_lowercase().as_str(), "y" | "yes")
}

/// Print a summary table of the GalaRP runs performed in `--all` mode.
fn print_summary(results: &[(&Plan, Outcome)], dry_run: bool, total: Duration) {
    let name = |p: &Plan| p.path.file_name()
        .map_or_else(|| p.path.display().to_string(), |n| n.to_string_lossy().into_owned());
    let width = results.iter().map(|(p, _)| name(p).len()).max().unwrap_or(0).max(6);
    let rule = "-".repeat(width + 22);

    println!("\n{:<width$}  {:<8}  {:>10}", "config", "status", "time");
    println!("{rule}");
    for (plan, outcome) in results {
        let (status, time) = match outcome {
            Outcome::Ran(d) => (if dry_run { "checked" } else { "ok" }, format!("{:.1} s", d.as_secs_f64())),
            Outcome::Skipped => ("skipped", "-".into()),
            Outcome::Failed(_) => ("FAILED", "-".into()),
        };
        println!("{:<width$}  {:<8}  {:>10}", name(plan), status, time);
    }
    println!("{rule}");

    let count = |f: fn(&Outcome) -> bool| results.iter().filter(|(_, o)| f(o)).count();
    println!("{} ok, {} skipped, {} failed in {:.1} s",
        count(|o| matches!(o, Outcome::Ran(_))),
        count(|o| matches!(o, Outcome::Skipped)),
        count(|o| matches!(o, Outcome::Failed(_))),
        total.as_secs_f64());

    // Repeat failure reasons here so they aren't lost in the scrollback.
    for (plan, outcome) in results {
        if let Outcome::Failed(e) = outcome {
            println!("  {}: {e}", name(plan));
        }
    }
}
