
use std::thread;
use std::time::{Duration, Instant};

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::sync::Arc;

use galarp_io::config::Config;
use galarp_run::{Run, GALAXY_COLUMNS};

pub enum Msg {
    Started { summary: Vec<String>, chunks: u64, n_clouds: usize},
    Step(Box<Step>),
    Finished(Result<bool, String>),

}

pub struct Step {
    pub t_myr: f64,
    pub row: Vec<f64>,

    pub snapshot: Option<u64>,

    pub wind_dir: Option<[f64; 3]>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Status {
    Running,
    Paused,
    Done,
    Stopped,
    Failed(String),
}

pub struct RunData {
    pub name: String,
    pub cfg: Config,
    pub status: Status,
    pub summary: Vec<String>,
    pub chunks: u64,
    pub n_clouds: usize,
    pub steps: Vec<Step>,
    pub started: Instant,
    pub elapsed: Duration,
    rx: Option<Receiver<Msg>>,
    stop: Arc<AtomicBool>,
}


impl RunData {
    pub fn start(name: String, cfg: Config, text: String) -> Self {
        let (tx, rx) = mpsc::channel();

        let stop = Arc::new(AtomicBool::new(false));
        let stop_worker = Arc::clone(&stop);
        // Clone the config object to allow the worker to use it
        let worker_cfg = cfg.clone();
        
        thread::spawn(move || {
            let result = (|| -> Result<bool, String> {
                let mut run = Run::new(worker_cfg, text).map_err(|e| e.to_string())?;
                let chunks = run.chunks();
                let n_clouds = run.clouds.galaxy.iter().filter(|&&g| g == 0).count();

                let _ = tx.send(Msg::Started { summary: run.summary.clone(), chunks, n_clouds });
                run.execute(|p| {
                    let step = Step { 
                        t_myr: p.t_myr, 
                        row: p.row.to_vec(), 
                        snapshot: p.snapshot, 
                        wind_dir: p.wind_dir[0] 
                    };
                    // A failed send means the UI dropped the run: stop.
                    if tx.send(Msg::Step(Box::new(step))).is_err() || stop_worker.load(Ordering::Relaxed) {
                        std::ops::ControlFlow::Break(())
                    } else {
                        std::ops::ControlFlow::Continue(())
                    }
                })
                .map_err(|e| e.to_string())
            })();
            let _ = tx.send(Msg::Finished(result));
        });

        Self {
            name,
            cfg,
            status: Status::Running,
            summary: Vec::new(),
            chunks: 0,
            n_clouds: 0,
            steps: Vec::new(),
            started: Instant::now(),
            elapsed: Duration::ZERO,
            rx: Some(rx),
            stop
        }
    }

    pub fn stop(&self) {
        self.stop.store(true, Ordering::Relaxed);
    }


    pub fn running(&self) -> bool {
        self.status == Status::Running
    }

    /// Take what the worker has sent since the last call.
    pub fn poll(&mut self) {
        let Some(rx) = &self.rx else { return };

        loop {
            match rx.try_recv() {
                Ok(Msg::Started { summary, chunks, n_clouds }) => {
                    (self.summary, self.chunks, self.n_clouds) = (summary, chunks, n_clouds)
                },
                Ok(Msg::Step(step)) => self.steps.push(*step),
                Ok(Msg::Finished(result)) => {
                    self.status = match result {
                        Ok(true) => Status::Done,
                        Ok(false) => Status::Stopped,
                        Err(e) => Status::Failed(e),
                    };
                    self.elapsed = self.started.elapsed();
                    self.rx = None;
                    return;
                }
                Err(TryRecvError::Empty) => {
                    self.elapsed = self.started.elapsed();
                    return;
                }
                Err(TryRecvError::Disconnected) => {
                    self.status = Status::Failed("run thread stopped unexpectedly".into());
                    self.rx = None;
                    return;
                }
            }
        }
    }

    /// Log column `name` of galaxy 0 at step `i`
    pub fn column(&self, i: usize, name: &str) -> f64 {
        let j = GALAXY_COLUMNS.iter().position(|&c| c == name).expect("a log column");
        self.steps[i].row[1 + j]
    }

}


impl Drop for RunData {
    fn drop(&mut self) {
        self.stop();
    }
}


// S in [0, 1] as a byte.
fn quantize(s: f64) -> u8 {
    (s.clamp(0.0, 1.0) * 255.0).round() as u8
}
