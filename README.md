
# GalaRP

**Idealized wind-tunnel simulations of ram-pressure stripping.**

Initially built on the Gala galaxy dynamics package ([Gala galaxy dynamics package](https://gala.adrian.pw/en/latest/)),  GalaRP models a galaxy's gas disk as test-particle clouds that move through analytic galaxy and cluster potentials. 

Each galaxy follows a precomputed orbit through the intracluster medium, and the clouds feel a ram-pressure wind with column-density shadowing. GalaRP is designed to be FAST, allowing for rapid traversal of the RPS parameter space.

## Installation
 
GalaRP needs Rust 1.85 or newer and the HDF5 C library.
 
```sh
# HDF5
brew install hdf5            # macOS
sudo apt install libhdf5-dev # Debian/Ubuntu
 
# The command-line tools
cargo install galarp-run galarp-tui
```

This verison of GalaRP takes in `.toml` config files (an example is in `configs/` ). and outputs a set of 
`hdf5` snaphot files that are set up to be readible by  ([`yt` ](https://yt-project.org/)) out of the box,  meaning you can 
immediately start running analysis.


## Building from Source

```sh
git clone https://github.com/HSouch/galarp
cd galarp
cargo build --release
cargo test
```

## Citing
A code paper is currently in development. Until then, please cite the first use of GalaRP in
Souchereau, et al (2026): DOI: 10.3847/1538-4357/ae976a 


## Licensing and Contributing

This code is available as open source. Please feel free to modify it if you see fit. Furthermore,
I welcome any suggested changes, modifications, and features. 