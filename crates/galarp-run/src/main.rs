//! Run GalaRP stripping sim from a TOML config.
//! 

use std::ops::ControlFlow;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use bigwig_io::config::Config;
use bigwig_io::snapshot::output_path;
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
}


fn main() -> ExitCode {
    let args = Args::parse();

    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> Result<(), Error> {
    let cfg = Config::from_file(&args.config)?;
    let config_text = std::fs::read_to_string(&args.config)?;

    let out = cfg.run.output_dir.clone();

    if !args.dry_run && output_path(&out, 0).exists() && !args.overwrite {
        return Err(
            format!("{} already has snapshots. Pass --overwrite or change run.output_dir", 
            out.display()).into() );
    }
    let t0 = Instant::now();
    println!("config: {}", args.config.display());

    let mut run = Run::new(cfg, config_text)?;

    for line in &run.summary {
        println!("{line}");
    }

    if args.dry_run {
        return Ok(());
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

    Ok(())

}