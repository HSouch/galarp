//! Tests for the whole config: parsing, defaults, validation and conversion.

use super::*;
use galarp::potentials::{Hernquist, Potential};
use galarp::units::{self, MYR, UNIT_DENSITY_CGS};

const MINIMAL: &str = r#"
    [run]
    output_dir = "out"
    t_end_myr = 100

    [box]
    side_kpc = 80

    [cluster.gas]
    model = "uniform"
    density_g_cm3 = 1e-27
    temperature_k = 5e7

    [[galaxy]]
    name = "parent"

    [[galaxy.potential]]
    type = "hernquist"
    mass_msun = 1e11
    a_kpc = 1.0

    [galaxy.disk]
    mass_msun = 5e9
    scale_length_kpc = 4
    r_max_kpc = 20
    sound_speed_km_s = 20
    n_particles = 5000

    [galaxy.orbit]
    mode = "constant_velocity"
    position_kpc = [0, 0, 0]
    velocity_km_s = [0, 0, 1000]
"#;

/// A second galaxy appended to MINIMAL.
const COMPANION: &str = r#"
    [[galaxy]]
    name = "companion"

    [[galaxy.potential]]
    type = "hernquist"
    mass_msun = 2e10
    a_kpc = 0.5

    [galaxy.disk]
    mass_msun = 1e9
    scale_length_kpc = 2
    r_max_kpc = 8
    sound_speed_km_s = 15
    n_particles = 300

    [galaxy.orbit]
    mode = "constant_velocity"
    position_kpc = [20, 0, 0]
    velocity_km_s = [0, 200, 1000]
"#;

fn parse(text: &str) -> Result<Config, ConfigError> {
    Config::from_str_in(text, Path::new("."))
}

#[test]
fn minimal_config_parses_with_defaults() {
    let c = parse(MINIMAL).unwrap();
    assert_eq!(c.run.dt_max_myr, 2.0);
    assert_eq!(c.snapshot_every(), 5);
    assert_eq!(c.hydro.gamma, units::GAMMA);
    assert_eq!(c.rk4_steps_per_dt_max().unwrap(), 4);
    assert!(c.cluster_potential().is_none());
    assert!(c.ram_pressure().unwrap().is_none(), "ram pressure is off by default");
    let t = c.orbit_table(0).unwrap();
    assert!(t.t_end() >= c.t_end());
    assert_eq!(c.galaxy_name(0), "parent");
    assert_eq!(c.disk(0).unwrap().1, 5000);
    match c.cluster_gas().unwrap() {
        ClusterGas::Uniform { rho, .. } => assert!((rho * UNIT_DENSITY_CGS / 1e-27 - 1.0).abs() < 1e-12),
        _ => panic!("expected uniform gas"),
    }
}

#[test]
fn unknown_keys_and_bad_values_are_errors() {
    let typo = MINIMAL.replace("side_kpc = 80", "sid_kpc = 80");
    let e = parse(&typo).unwrap_err();
    assert!(e.0.contains("sid_kpc") || e.0.contains("side_kpc"), "{e}");
    let unknown_type = MINIMAL.replace("type = \"hernquist\"", "type = \"plummer\"");
    assert!(parse(&unknown_type).is_err());
    let bad = MINIMAL.replace("t_end_myr = 100", "t_end_myr = -1");
    assert!(parse(&bad).unwrap_err().0.contains("t_end_myr"));
    let misaligned = MINIMAL.replace("[run]", "[run]\nsnapshot_interval_myr = 3");
    assert!(parse(&misaligned).unwrap_err().0.contains("multiple"));
    let beta_without_potential = MINIMAL.replace(
        "model = \"uniform\"\n    density_g_cm3 = 1e-27\n    temperature_k = 5e7",
        "model = \"beta\"\n    rho0_g_cm3 = 1e-26\n    r_c_kpc = 100\n    beta = 0.6",
    );
    assert!(parse(&beta_without_potential).unwrap_err().0.contains("cluster.potential"));
    let no_count = MINIMAL.replace("n_particles = 5000", "");
    assert!(parse(&no_count).unwrap_err().0.contains("n_particles"));
}

/// SPH-only sections are not part of this schema.
#[test]
fn full_sph_sections_are_rejected() {
    for extra in ["[self_gravity]\nenabled = true\n", "[cooling]\nenabled = true\n", "[glass]\nn_particles = 100\n", "[feedback]\nenabled = true\n"] {
        assert!(parse(&format!("{MINIMAL}\n{extra}")).is_err(), "{extra}");
    }
    assert!(parse(&MINIMAL.replace("side_kpc = 80", "side_kpc = 80\n    recenter_fraction = 0.25")).is_err());
    assert!(parse(&MINIMAL.replace("[box]", "[hydro]\n    n_ngb = 64\n\n    [box]")).is_err());
}

#[test]
fn time_dependent_parameters_and_unit_conversion() {
    let text = MINIMAL.replace(
        "mass_msun = 1e11\n    a_kpc = 1.0",
        "mass_msun = { t_myr = [0, 100], value = [1e11, 5e10] }\n    a_kpc = 1.0",
    );
    let c = parse(&text).unwrap();
    let p = c.galaxy_potential(0);
    let x = [50.0, 0.0, 0.0];
    let a0 = p.acceleration(&x, 0.0)[0];
    let a1 = p.acceleration(&x, 100.0 * MYR)[0];
    assert!((a1 / a0 - 0.5).abs() < 1e-12);
    let expected = Hernquist { mass: 10.0.into(), a: 1.0.into() }.acceleration(&x, 0.0)[0];
    assert!((a0 / expected - 1.0).abs() < 1e-12);
}

#[test]
fn ram_pressure_and_galarp_sections() {
    let text = format!("{MINIMAL}\n[ram_pressure]\nenabled = true\nsurface_density_msun_pc2 = 20\n[ram_pressure.shadow]\nenabled = true\n[galarp]\ndt_myr = 0.25\n");
    let c = parse(&text).unwrap();
    let rp = c.ram_pressure().unwrap().unwrap();
    assert!(rp.shadow.is_some());
    assert_eq!(c.rk4_steps_per_dt_max().unwrap(), 8);
    assert!(parse(&text.replace("dt_myr = 0.25", "dt_myr = 0.3")).unwrap().rk4_steps_per_dt_max().is_err());
    assert!(parse(&text.replace("surface_density_msun_pc2 = 20", "surface_density_msun_pc2 = -1")).is_err());
    let none = text.replace("model = \"uniform\"\n    density_g_cm3 = 1e-27\n    temperature_k = 5e7", "model = \"none\"");
    assert!(parse(&none).unwrap_err().0.contains("ram_pressure"));
}

#[test]
fn disk_options_from_the_disk_section() {
    let c = parse(MINIMAL).unwrap();
    let o = c.disk_options(0);
    assert!(!o.self_gravity && o.noise_fraction == 0.0);
    let sg = MINIMAL.replace("n_particles = 5000", "n_particles = 5000\n    include_self_gravity = true\n    noise_fraction = 0.1");
    let o2 = parse(&sg).unwrap().disk_options(0);
    assert!(o2.self_gravity);
    assert_eq!(o2.noise_fraction, 0.1);
}

#[test]
fn galaxy_orientation_from_a_normal_or_angles() {
    let with = |extra: &str| MINIMAL.replacen("[[galaxy.potential]]", &format!("{extra}\n[[galaxy.potential]]"), 1);
    let c = parse(MINIMAL).unwrap();
    assert!(c.galaxy_orientation(0).is_identity());
    let n = parse(&with("disk_normal = [0, 1, 1]")).unwrap().galaxy_orientation(0).normal();
    assert!((n[1] - 0.5f64.sqrt()).abs() < 1e-12 && (n[2] - 0.5f64.sqrt()).abs() < 1e-12);
    let n = parse(&with("inclination_deg = 90\nazimuth_deg = 90")).unwrap().galaxy_orientation(0).normal();
    assert!((n[1] - 1.0).abs() < 1e-12);
    assert!(parse(&with("disk_normal = [0, 0, 1]\ninclination_deg = 10")).is_err(), "both");
    assert!(parse(&with("disk_normal = [0, 0, 0]")).is_err(), "zero");
    assert!(parse(&with("azimuth_deg = 10")).is_err(), "azimuth alone");
    assert!(parse(&with("inclination_deg = 200")).is_err(), "out of range");
}

#[test]
fn no_icm_model_parses() {
    let text = MINIMAL.replace("model = \"uniform\"\n    density_g_cm3 = 1e-27\n    temperature_k = 5e7", "model = \"none\"");
    let c = parse(&text).unwrap();
    assert!(matches!(c.cluster_gas().unwrap(), ClusterGas::Uniform { rho, .. } if rho == 0.0));
}

#[test]
fn a_second_galaxy_has_its_own_cloud_count() {
    let c = parse(&format!("{MINIMAL}{COMPANION}")).unwrap();
    assert_eq!(c.galaxy.len(), 2);
    assert_eq!(c.galaxy_name(1), "companion");
    assert_eq!(c.disk(1).unwrap().1, 300, "clouds need not share a mass");
    let p0 = c.orbit_table(0).unwrap().at(0.0).unwrap().0;
    let p1 = c.orbit_table(1).unwrap().at(0.0).unwrap().0;
    assert_eq!((p0[0], p1[0]), (0.0, 20.0));
    // Without any galaxy the config is rejected.
    let none = MINIMAL.split("[[galaxy]]").next().unwrap();
    assert!(parse(none).is_err());
}

#[test]
fn integrated_galaxies_pull_on_each_other() {
    // Two equal Hernquist galaxies starting at rest 40 kpc apart, with no cluster: they
    // fall together symmetrically. Without mutual gravity nothing would pull them.
    let pair = |x: f64| {
        format!(
            "[[galaxy]]\nname = \"g{x}\"\n[[galaxy.potential]]\ntype = \"hernquist\"\nmass_msun = 1e11\na_kpc = 1.0\n[galaxy.orbit]\nmode = \"integrate\"\nposition_kpc = [{x}, 0, 0]\nvelocity_km_s = [0, 0, 0]\n"
        )
    };
    let base = MINIMAL.split("[[galaxy]]").next().unwrap().replace("t_end_myr = 100", "t_end_myr = 200");
    let c = parse(&format!("{base}{}{}", pair(20.0), pair(-20.0))).unwrap();
    let t = c.orbit_tables().unwrap();
    let (a, b) = (t[0].at(c.t_end()).unwrap().0, t[1].at(c.t_end()).unwrap().0);
    assert!(a[0] < 19.0 && (a[0] + b[0]).abs() < 1e-9, "{a:?} {b:?}");
    let off = format!("{base}[orbits]\nmutual_gravity = false\n{}{}", pair(20.0), pair(-20.0));
    assert!(parse(&off).unwrap_err().0.contains("cluster.potential"));
}

#[test]
fn integrated_orbit_needs_and_uses_the_cluster_potential() {
    let text = MINIMAL
        .replace("mode = \"constant_velocity\"", "mode = \"integrate\"")
        .replace("[cluster.gas]", "[[cluster.potential]]\n    type = \"nfw_virial\"\n    m_vir_msun = 1e15\n    r_vir_kpc = 2000\n    concentration = 5\n\n    [cluster.gas]")
        .replace("position_kpc = [0, 0, 0]", "position_kpc = [1000, 0, 0]")
        .replace("velocity_km_s = [0, 0, 1000]", "velocity_km_s = [0, 300, 0]");
    let c = parse(&text).unwrap();
    let t = c.orbit_table(0).unwrap();
    let (x, _) = t.at(c.t_end()).unwrap();
    assert!(x[0] < 1000.0, "galaxy should fall inward");
}

#[test]
fn typo_inside_a_potential_term_is_reported_with_its_location() {
    let text = MINIMAL.replace("a_kpc = 1.0", "a_kcp = 1.0");
    let e = parse(&text).unwrap_err();
    assert!(e.0.contains("a_kcp"), "{e}");
}
