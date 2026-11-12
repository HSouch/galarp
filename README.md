# GalaRP

**Idealized wind-tunnel simulations of ram-pressure stripping.**

Initially built on the [Gala galaxy dynamics package](https://gala.adrian.pw/en/latest/), GalaRP models a galaxy's gas disk as test-particle clouds that move through analytic galaxy and cluster potentials.

Each galaxy follows a precomputed orbit through the intracluster medium, and the clouds feel a ram-pressure wind with column-density shadowing. GalaRP is designed to be FAST, allowing for rapid traversal of the RPS parameter space.

## Crates

| Crate | Function |
|---|---|
| [`galarp`](https://crates.io/crates/galarp) | Core library: potentials, cluster profiles, orbits, initial conditions, and the ram-pressure integrator |
| [`galarp-io`](https://crates.io/crates/galarp-io) | TOML config parsing, HDF5 snapshot output, and CSV logging |
| [`galarp-run`](https://crates.io/crates/galarp-run) | Run simulations from a config file in the terminal |
| [`galarp-tui`](https://crates.io/crates/galarp-tui) | Terminal UI to make, edit, run, and watch runs |

## Installation

GalaRP needs Rust 1.88 or newer and the HDF5 C library.

```sh
# HDF5
brew install hdf5            # macOS
sudo apt install libhdf5-dev # Debian/Ubuntu

# The command-line tools
cargo install galarp-run galarp-tui
```

To use the library in your own code:

```sh
cargo add galarp galarp-io
```

## Usage

This version of GalaRP takes in `.toml` config files (examples are in [`configs/`](https://github.com/HSouch/galarp/tree/master/configs)) and outputs a set of
`hdf5` snapshot files that are set up to be readable by [`yt`](https://yt-project.org/) out of the box, meaning you can
immediately start running analysis.

```sh
galarp-run configs/face_on.toml
```

## Building from Source

```sh
git clone https://github.com/HSouch/galarp
cd galarp
cargo build --release
cargo test
```

## Citing

A code paper is currently in development. Until then, please cite the first use of GalaRP in
Souchereau, et al (2026): DOI: [10.3847/1538-4357/ae976a](https://doi.org/10.3847/1538-4357/ae976a)

[ADS Link](https://ui.adsabs.harvard.edu/abs/2026ApJ..1008..182S/abstract)

## License

Licensed under the MIT license.

- MIT license ([LICENSE-MIT](LICENSE-MIT) or <https://opensource.org/licenses/MIT>)


## Contributing

I welcome any suggested changes, modifications, and features. Feel free to submit pull requests,
or contact me about ideas.
