//! HDF5 snapshot saving.
//! Should be readable by yt by default, as the setup is similar to GADGET-4.
//! 
//! Layout:
//!     Header
//!     Parameters
//!     Config
//!     Clouds (PartType1) particle objects
//! 

use std::path::{Path, PathBuf};

use galarp::analytic::Clouds;
use galarp::units;
use hdf5::types::VarLenUnicode;
use hdf5::{File, Group, H5Type};

pub use hdf5::Result;

/// Number of GADGET particle types.
const N_TYPES: usize = 6;


#[derive(Clone, Debug, PartialEq)]
pub struct SnapshotHeader {
    pub time: f64,
    /// Side of the cubic box (code length units).
    pub box_size: f64,
    /// Position of the box-frame origin in cluster coordinates.
    pub box_origin: [f64; 3],
    /// Center position and velocity of every galaxy in cluster coordinates; index 0 is the
    /// parent galaxy.
    pub galaxy_pos: Vec<[f64; 3]>,
    pub galaxy_vel: Vec<[f64; 3]>,
    /// Each galaxy's disk spin axis in the cluster frame (empty: all [0, 0, 1]).
    pub galaxy_normal: Vec<[f64; 3]>,
    /// Adiabatic index, and the mean molecular weight used for the Temperature field.
    pub gamma: f64,
    pub mean_molecular_weight: f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct WriteOptions {
    /// Apply shuffle + gzip to particle datasets.
    pub compress: bool,
}

/// Gen new attribute
fn attr<T: H5Type>(g: &Group, name: &str, value: T) -> Result<()> {
    g.new_attr::<T>().create(name)?.write_scalar(&value)
}

/// Gen new array attribute
fn attr_array<T: H5Type>(g: &Group, name: &str, values: &[T]) -> Result<()> {
    g.new_attr_builder().with_data(values).create(name)?;
    Ok(())
}

/// Gen new string attribute
fn attr_string(g: &Group, name: &str, value: &str) -> Result<()> {
    let v: VarLenUnicode = value.parse().map_err(|e| hdf5::Error::from(format!("{e:?}")))?;
    g.new_attr::<VarLenUnicode>().create(name)?.write_scalar(&v)
}

fn read_attr<T: H5Type>(g: &Group, name: &str) -> Result<T> {
    g.attr(name)?.read_scalar()
}

/// 1D Dataset
fn dataset_1d<T: H5Type>(g: &Group, name: &str, data: &[T], opts: WriteOptions) -> Result<()> {
    let b = g.new_dataset_builder().with_data(data);
    if opts.compress && !data.is_empty() {
        b.chunk(data.len().min(1 << 16)).shuffle().deflate(4).create(name)?;
    } else {
        b.create(name)?;
    }
    Ok(())
}


fn dataset_3d(g: &Group, name: &str, data: &[[f64; 3]], opts: WriteOptions) -> Result<()> {
    let b = g.new_dataset::<f64>().shape((data.len(), 3));
    let ds = if opts.compress && !data.is_empty() { 
        b.chunk((data.len().min(1 << 16), 3))
        .shuffle().deflate(4).create(name)? 

    } else { b.create(name)? };
    ds.write_raw(data.as_flattened())
}

/// Name of output `number`, Enzo style: `DD0012`. Each output is a directory of that name
/// holding one file, `DD0012/DD0012.hdf5`.
pub fn output_name(number: u64) -> String {
    format!("DD{number:04}")
}

/// Path of output `number` in run directory `run_dir`: `run_dir/DD0012/DD0012.hdf5`.
pub fn output_path(run_dir: &Path, number: u64) -> PathBuf {
    let name = output_name(number);
    run_dir.join(&name).join(format!("{name}.hdf5"))
}

/// The snapshot file for a restart argument: an output file, or an output directory
/// (`DD0012/`), in which case the file inside it.
pub fn resolve_output(path: &Path) -> PathBuf {
    if path.is_dir() {
        if let Some(name) = path.file_name() {
            return path.join(format!("{}.hdf5", name.to_string_lossy()));
        }
    }
    path.to_path_buf()
}

/// The run directory a snapshot file belongs to: two levels up for `DD0012/DD0012.hdf5`,
/// one for a file directly in the run directory (the older `snap_012.hdf5` layout).
pub fn run_dir_of(snapshot: &Path) -> PathBuf {
    let parent = snapshot.parent().unwrap_or(Path::new("."));
    let in_own_dir = snapshot.file_stem().zip(parent.file_name()).is_some_and(|(stem, dir)| stem == dir && stem.to_string_lossy().starts_with("DD"));
    let dir = if in_own_dir { parent.parent().unwrap_or(Path::new(".")) } else { parent };
    if dir.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        dir.to_path_buf()
    }
}

/// The `Header`, `Parameters` and `Config` groups, with `num_part` particles of each type.
fn write_header_groups(file: &File, header: &SnapshotHeader, num_part: &[u32; N_TYPES]) -> Result<()> {
    let h = file.create_group("Header")?;
    attr_array(&h, "NumPart_ThisFile", num_part)?;
    attr_array(&h, "NumPart_Total", num_part)?;
    attr_array(&h, "NumPart_Total_HighWord", &[0u32; N_TYPES])?;
    attr_array(&h, "MassTable", &[0.0f64; N_TYPES])?;
    attr(&h, "Time", header.time)?;
    attr(&h, "Redshift", 0.0f64)?;
    attr(&h, "BoxSize", header.box_size)?;
    attr(&h, "NumFilesPerSnapshot", 1i32)?;
    attr(&h, "Flag_DoublePrecision", 1i32)?;
    attr_array(&h, "BoxOrigin", &header.box_origin)?;
    let n_gal = header.galaxy_pos.len();
    assert_eq!(n_gal, header.galaxy_vel.len());
    attr(&h, "NumGalaxies", n_gal as i32)?;
    let parent_pos = header.galaxy_pos.first().copied().unwrap_or([0.0; 3]);
    let parent_vel = header.galaxy_vel.first().copied().unwrap_or([0.0; 3]);
    attr_array(&h, "ParentGalaxyPosition", &parent_pos)?;
    attr_array(&h, "ParentGalaxyVelocity", &parent_vel)?;
    if n_gal > 0 {
        attr_array(&h, "GalaxyPositions", header.galaxy_pos.as_flattened())?;
        attr_array(&h, "GalaxyVelocities", header.galaxy_vel.as_flattened())?;
        let normals: Vec<[f64; 3]> = (0..n_gal).map(|k| header.galaxy_normal.get(k).copied().unwrap_or([0.0, 0.0, 1.0])).collect();
        attr_array(&h, "GalaxyDiskNormals", normals.as_flattened())?;
        attr_array(&h, "ParentDiskNormal", &normals[0])?;
    }

    let p = file.create_group("Parameters")?;
    attr(&p, "UnitLength_in_cm", units::UNIT_LENGTH_CM)?;
    attr(&p, "UnitMass_in_g", units::UNIT_MASS_G)?;
    attr(&p, "UnitVelocity_in_cm_per_s", units::UNIT_VELOCITY_CM_S)?;
    attr(&p, "Gamma", header.gamma)?;
    attr(&p, "MeanMolecularWeight", header.mean_molecular_weight)?;

    // Must not contain VORONOI or AMR, or yt would read the file as Arepo.
    let c = file.create_group("Config")?;
    attr(&c, "BIGWIG_VERSION_MAJOR", env!("CARGO_PKG_VERSION_MAJOR").parse::<i32>().unwrap())?;
    Ok(())
}

pub fn write_cloud_snapshot(
    path: impl AsRef<Path>,
    header: &SnapshotHeader,
    clouds: &Clouds,
    galaxy_coords: &[[f64; 3]],
    config_text: &str,
    opts: WriteOptions,
) -> Result<()> {
    assert_eq!(galaxy_coords.len(), clouds.len());
    let file = File::create(path)?;
    let mut num_part = [0u32; N_TYPES];
    num_part[1] = u32::try_from(clouds.len()).expect("more than 2^32 clouds");
    write_header_groups(&file, header, &num_part)?;
    let eps = 1e-9 * header.box_size;
    let o = header.box_origin;
    let pos_out: Vec<[f64; 3]> =
        clouds.pos.iter().map(|x| [x[0] - o[0], x[1] - o[1], x[2] - o[2]].map(|c| c.clamp(0.0, header.box_size - eps))).collect();
    let g = file.create_group("PartType1")?;
    dataset_3d(&g, "Coordinates", &pos_out, opts)?;
    dataset_3d(&g, "ClusterCoordinates", &clouds.pos, opts)?;
    dataset_3d(&g, "GalaxyCoordinates", galaxy_coords, opts)?;
    dataset_3d(&g, "Velocities", &clouds.vel, opts)?;
    dataset_1d(&g, "ParticleIDs", &clouds.id, opts)?;
    dataset_1d(&g, "Masses", &clouds.mass, opts)?;
    let origin: Vec<i32> = clouds.galaxy.iter().map(|&k| k as i32).collect();
    dataset_1d(&g, "Origin", &origin, opts)?;
    dataset_1d(&g, "Shadowing", &clouds.shadow, opts)?;
    dataset_3d(&g, "RamPressureAcceleration", &clouds.a_ram, opts)?;
    let r = file.create_group("Restricted")?;
    attr_string(&r, "ConfigToml", config_text)?;
    Ok(())
}



/// The clouds of a  `galaRP` snapshot, as `write_cloud_snapshot` wrote them (the
/// parts a viewer needs).
#[derive(Clone, Debug, PartialEq)]
pub struct CloudSnapshot {
    /// Time of the output (code units).
    pub time: f64,
    /// Each cloud's galaxy (`Origin`), its offset from that galaxy's center in the galaxy's
    /// frame (`GalaxyCoordinates`), and its shadow factor (`Shadowing`).
    pub galaxy: Vec<u16>,
    pub galaxy_coords: Vec<[f64; 3]>,
    pub shadow: Vec<f64>,
}

/// Reads a `bigwig-restricted` snapshot (a `DDNNNN` directory or its file).
pub fn read_cloud_snapshot(path: impl AsRef<Path>) -> Result<CloudSnapshot> {
    let file = File::open(resolve_output(path.as_ref()))?;
    let g = file.group("PartType1")?;
    let flat: Vec<f64> = g.dataset("GalaxyCoordinates")?.read_raw()?;
    let origin: Vec<i32> = g.dataset("Origin")?.read_raw()?;
    Ok(CloudSnapshot {
        time: read_attr(&file.group("Header")?, "Time")?,
        galaxy: origin.into_iter().map(|k| k as u16).collect(),
        galaxy_coords: flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect(),
        shadow: g.dataset("Shadowing")?.read_raw()?,
    })
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn outputs_are_named_and_found_enzo_style() {
        let p = output_path(Path::new("runs/x"), 12);
        assert_eq!(p, Path::new("runs/x/DD0012/DD0012.hdf5"));
        assert_eq!(run_dir_of(&p), Path::new("runs/x"));
        assert_eq!(run_dir_of(Path::new("DD0003/DD0003.hdf5")), Path::new("."));
        // The older layout: snapshots directly in the run directory.
        assert_eq!(run_dir_of(Path::new("runs/x/snap_012.hdf5")), Path::new("runs/x"));
        assert_eq!(run_dir_of(Path::new("snap_012.hdf5")), Path::new("."));
        assert_eq!(output_name(12345), "DD12345");
    }
}
