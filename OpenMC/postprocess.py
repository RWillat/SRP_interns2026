import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openmc

# Shared postprocessing for the HPMR_pinCell (TRISO), HPMR_pinCell-UN (solid
# UN + SiC clad), and HPMR_INFINIT (UN/YH2 heat pipe) unit-cell models: per
# case, reaction-rate/heating tables, flux spectra, flux/power-density mesh
# heatmaps, and a specific-power (fission power per loaded HM mass) metric;
# then cross-case spectrum and specific-power comparisons.
#
# OpenMC eigenvalue tallies are per source neutron. Each case is normalized
# so its total recoverable heat equals P0_WATTS (arbitrary; results scale
# linearly with it), which makes specific power intensive (W per gram HM,
# independent of the modeled domain size). A fair cross-case comparison
# still assumes both cells run at the same reference power P0.

HERE = os.path.dirname(os.path.abspath(__file__))

P0_WATTS = 1.0  # reference total recoverable power deposited per modeled cell

# Per-case config: fuel material name and the heavy-metal nuclides that make
# up the "loaded uranium / heavy metal" denominator of the specific-power
# metric. get_mass_density() must be summed one nuclide at a time (passing a
# list returns 0.0 in openmc 0.15.x).
CASES = {
    'HPMR_pinCell': {'fuel': 'fuel', 'hm_nuclides': ['U235', 'U238']},
    'HPMR_pinCell-UN': {'fuel': 'UN', 'hm_nuclides': ['U235', 'U238']},
    'HPMR_INFINIT': {'fuel': 'UN', 'hm_nuclides': ['U235', 'U238']},
}


def mat_labels(material_filter, material_names):
    return [material_names[int(m)] for m in material_filter.bins]


def load_case(case_dir):
    statepoint_path = sorted(glob.glob(os.path.join(case_dir, 'statepoint.*.h5')))[-1]
    sp = openmc.StatePoint(statepoint_path)
    material_names = {m.id: m.name for m in sp.summary.materials}
    print(f'Loaded {statepoint_path}')
    print(f'k-eff: {sp.keff}\n')
    return sp, material_names


def reaction_rate_tables(case_dir, sp, material_names):
    absorption = sp.get_tally(name='absorption-by-material')
    mat_filter = absorption.find_filter(openmc.MaterialFilter)
    mean = absorption.get_reshaped_data(value='mean').reshape(len(mat_filter.bins), -1)[:, 0]
    std = absorption.get_reshaped_data(value='std_dev').reshape(len(mat_filter.bins), -1)[:, 0]
    df_abs = pd.DataFrame({
        'material': mat_labels(mat_filter, material_names),
        'absorption rate': mean,
        'std. dev.': std,
    })
    print('=== Absorption rate by material ===')
    print(df_abs.to_string(index=False))
    print()
    df_abs.to_csv(os.path.join(case_dir, 'absorption_by_material.csv'), index=False)

    fission = sp.get_tally(name='fission-in-fuel')
    mean = fission.get_reshaped_data(value='mean').reshape(-1)
    std = fission.get_reshaped_data(value='std_dev').reshape(-1)
    df_fis = pd.DataFrame({'score': fission.scores, 'mean': mean, 'std. dev.': std})
    print('=== Fission / nu-fission / kappa-fission rate in fuel ===')
    print(df_fis.to_string(index=False))
    print()
    df_fis.to_csv(os.path.join(case_dir, 'fission_in_fuel.csv'), index=False)


def heating_by_material(case_dir, sp, material_names):
    """Local heating rate per material, with each material's share of the
    total recoverable heat deposited in the cell."""
    heating = sp.get_tally(name='heating-by-material')
    mat_filter = heating.find_filter(openmc.MaterialFilter)
    mean = heating.get_reshaped_data(value='mean').reshape(len(mat_filter.bins), -1)[:, 0]
    std = heating.get_reshaped_data(value='std_dev').reshape(len(mat_filter.bins), -1)[:, 0]
    total = mean.sum()
    df = pd.DataFrame({
        'material': mat_labels(mat_filter, material_names),
        'heating rate [eV/src]': mean,
        'std. dev.': std,
        'fraction of total': mean / total,
    }).sort_values('heating rate [eV/src]', ascending=False)
    print('=== Local heating rate by material (heating-local) ===')
    print(df.to_string(index=False))
    print(f'total heating: {total:.6e} eV/src\n')
    df.to_csv(os.path.join(case_dir, 'heating_by_material.csv'), index=False)
    return total


def specific_power(case_dir, sp, material_names, cfg, total_heating):
    """Fission power per unit loaded heavy-metal mass.

    Normalizes so total deposited power in the cell = P0_WATTS, attributes the
    recoverable fission energy (kappa-fission in fuel) as the fuel fission
    power, and divides by the heavy-metal mass loaded in the modeled cell.
    """
    fuel_name = cfg['fuel']
    fuel_mat = next(m for m in sp.summary.materials if m.name == fuel_name)
    if fuel_mat.volume is None:
        raise ValueError(f'{fuel_name} volume is not set in the statepoint; '
                         'set mats[...].volume in the model before running.')

    # HM mass density: sum single-nuclide get_mass_density calls (a list arg
    # returns 0.0 in openmc 0.15.x).
    hm_density = sum(fuel_mat.get_mass_density(nuc) for nuc in cfg['hm_nuclides'])
    hm_mass = hm_density * fuel_mat.volume  # g

    fission = sp.get_tally(name='fission-in-fuel')
    kappa = fission.get_values(scores=['kappa-fission']).ravel()[0]  # eV/src

    # Watts per (eV/src) so that the whole cell deposits P0_WATTS.
    norm = P0_WATTS / total_heating
    fuel_fission_power = kappa * norm  # W
    specific_power_w_per_kg = fuel_fission_power / (hm_mass / 1000)  # W/kgHM
    power_density_w_per_cc = fuel_fission_power / fuel_mat.volume  # W/cc fuel

    result = {
        'case': os.path.basename(case_dir),
        'fuel material': fuel_name,
        'HM density [g/cc]': hm_density,
        'fuel volume [cc]': fuel_mat.volume,
        'HM mass [g]': hm_mass,
        'kappa-fission [eV/src]': kappa,
        'fuel fission fraction of heat': kappa / total_heating,
        f'specific power [W/kgHM @ {P0_WATTS:g} W cell]': specific_power_w_per_kg,
        f'power density [W/cc @ {P0_WATTS:g} W cell]': power_density_w_per_cc,
    }
    print('=== Specific fission power (per loaded heavy metal) ===')
    for k, v in result.items():
        if k == 'case' or k == 'fuel material':
            continue
        print(f'  {k:45s}: {v:.6e}')
    print()
    return result


def get_spectra(sp, material_names):
    """Return (labels, e_bins, e_mid, lethargy_width, flux[n_mat, n_e])."""
    spectrum = sp.get_tally(name='flux-spectrum-by-region')
    mat_filter = spectrum.find_filter(openmc.MaterialFilter)
    energy_filter = spectrum.find_filter(openmc.EnergyFilter)
    labels = mat_labels(mat_filter, material_names)
    e_bins = np.array(energy_filter.bins)  # shape (n_groups, 2), eV
    n_mat, n_e = len(labels), len(e_bins)

    flux = spectrum.get_reshaped_data(value='mean').reshape(n_mat, n_e)

    e_lo, e_hi = e_bins[:, 0], e_bins[:, 1]
    e_mid = np.sqrt(np.clip(e_lo, 1e-5, None) * e_hi)  # log-midpoint
    lethargy_width = np.log(e_hi / np.clip(e_lo, 1e-5, None))
    return labels, e_bins, e_mid, lethargy_width, flux


def spectra_plot(case_dir, case_name, sp, material_names):
    labels, e_bins, e_mid, lethargy_width, flux = get_spectra(sp, material_names)
    flux_per_lethargy = flux / lethargy_width

    # Normalize each region's spectrum to unit area (in lethargy space) so
    # spectral *shape* can be compared across regions/cases independent of
    # each region's overall flux magnitude.
    area = np.sum(flux_per_lethargy * lethargy_width, axis=1, keepdims=True)
    flux_normalized = flux_per_lethargy / area

    for normalized, suffix, ylabel in [
        (False, '', r'Flux per unit lethargy [#/cm$^2\cdot u$]'),
        (True, '_normalized', 'Flux per Lethargy [Normalized]'),
    ]:
        data = flux_normalized if normalized else flux_per_lethargy
        fig, ax = plt.subplots(figsize=(8, 6))
        for i, name in enumerate(labels):
            ax.step(e_mid, data[i, :], where='mid', label=name)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel(ylabel)
        title = f'{case_name}: flux spectrum by region'
        ax.set_title(title + (' (normalized)' if normalized else ''))
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(case_dir, f'flux_spectra{suffix}.png'), dpi=150)
        plt.close(fig)

    df_spec = pd.DataFrame(flux_per_lethargy, index=labels,
                            columns=[f'{lo:.4g}-{hi:.4g} eV' for lo, hi in e_bins])
    df_spec.to_csv(os.path.join(case_dir, 'flux_spectrum_by_region.csv'))

    return labels, e_mid, flux_normalized


def mesh_heatmap(case_dir, case_name, sp, tally_name, out_name, title_quantity,
                  cbar_label, norm_factor=1.0):
    """Plot the x-y density map of a mesh tally's first (only) axial layer.

    OpenMC mesh-tally means are volume-integrated (extensive) per cell, not
    densities: summing 'heating-mesh' over all cells reproduces the
    non-meshed 'heating-by-material' total. Each cell must be divided by its
    own volume before plotting, otherwise the map's units mean nothing.
    """
    mesh_tally = sp.get_tally(name=tally_name)
    mesh_filter = mesh_tally.find_filter(openmc.MeshFilter)
    mesh = mesh_filter.mesh
    nx, ny, nz = mesh.dimension
    data = mesh_tally.get_reshaped_data(value='mean').reshape(nx, ny, nz)
    volumes = mesh.volumes.reshape(nx, ny, nz)
    density = data / volumes * norm_factor

    fig, ax = plt.subplots(figsize=(7, 6))
    extent = (mesh.lower_left[0], mesh.upper_right[0],
              mesh.lower_left[1], mesh.upper_right[1])
    # Mask cells with zero score (outside the model boundary) so the color
    # scale reflects the actual in-cell contrast instead of being stretched
    # by the zero background.
    plot_data = np.ma.masked_equal(density[:, :, 0].T, 0.0)
    cmap = plt.get_cmap('inferno').copy()
    cmap.set_bad('lightgray')
    im = ax.imshow(plot_data, origin='lower', extent=extent, cmap=cmap)
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_title(f'{case_name}: spatial {title_quantity} distribution')
    fig.colorbar(im, ax=ax, label=cbar_label)
    fig.tight_layout()
    fig.savefig(os.path.join(case_dir, out_name), dpi=150)
    plt.close(fig)


## --- Run per-case postprocessing, then combined cross-case comparisons --- ##

case_spectra = {}
sp_results = []
for case_name, cfg in CASES.items():
    case_dir = os.path.join(HERE, case_name)
    if not glob.glob(os.path.join(case_dir, 'statepoint.*.h5')):
        print(f'No statepoint found for {case_name}, skipping.\n')
        continue
    print(f'########## {case_name} ##########')
    sp, material_names = load_case(case_dir)
    reaction_rate_tables(case_dir, sp, material_names)
    total_heating = heating_by_material(case_dir, sp, material_names)
    sp_results.append(specific_power(case_dir, sp, material_names, cfg, total_heating))
    labels, e_mid, flux_normalized = spectra_plot(case_dir, case_name, sp, material_names)
    mesh_heatmap(case_dir, case_name, sp, 'flux-mesh', 'flux_mesh_heatmap.png',
                 title_quantity='flux density',
                 cbar_label='Flux [cm$^{-2}$ per source neutron]')
    mesh_heatmap(case_dir, case_name, sp, 'heating-mesh', 'heating_mesh_heatmap.png',
                 title_quantity='power density',
                 cbar_label=f'Power density [W/cm$^3$ @ {P0_WATTS:g} W cell]',
                 norm_factor=P0_WATTS / total_heating)
    case_spectra[case_name] = (labels, e_mid, flux_normalized)
    print(f'Wrote per-case outputs for {case_name}\n')

# Fuel-region spectra only, so the differing material sets per case stay
# readable on one plot
if len(case_spectra) > 1:
    fig, ax = plt.subplots(figsize=(8, 6))
    for case_name, (labels, e_mid, flux_normalized) in case_spectra.items():
        fuel_label = CASES[case_name]['fuel']
        if fuel_label not in labels:
            continue
        i = labels.index(fuel_label)
        ax.step(e_mid, flux_normalized[i, :], where='mid',
                label=f'{case_name} ({fuel_label})')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('Normalized flux per unit lethargy [1/eV, unit area]')
    ax.set_title('Fuel spectrum comparison across cases (normalized)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'fuel_spectrum_comparison_normalized.png'), dpi=150)
    plt.close(fig)
    print('Wrote fuel_spectrum_comparison_normalized.png')

# Cross-case specific-power comparison (table + bar chart)
if sp_results:
    df_sp = pd.DataFrame(sp_results)
    df_sp.to_csv(os.path.join(HERE, 'specific_power_comparison.csv'), index=False)
    print('\n=== Specific fission power comparison across cases ===')
    print(df_sp.to_string(index=False))
    print(f'\n(Normalized to P0 = {P0_WATTS:g} W total recoverable power per '
          'modeled cell; specific power scales linearly with P0.)')

    sp_col = f'specific power [W/kgHM @ {P0_WATTS:g} W cell]'
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.get_cmap('tab10')(np.linspace(0, 1, len(df_sp)))
    ax.bar(df_sp['case'], df_sp[sp_col], color=colors)
    ax.set_ylabel(sp_col)
    ax.set_yscale('log')
    ax.set_title('Fission power per unit loaded heavy-metal mass')
    for x, y in zip(df_sp['case'], df_sp[sp_col]):
        ax.annotate(f'{y:.3e}', (x, y), ha='center', va='bottom')
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'specific_power_comparison.png'), dpi=150)
    plt.close(fig)
    print('Wrote specific_power_comparison.png')
