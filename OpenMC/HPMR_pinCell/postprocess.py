import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openmc

# ---------------------------------------------------------------------------
# Postprocessing for HPMR_pincell.py: reaction-rate tables, region spectra,
# and a mesh flux heatmap, read from the eigenvalue statepoint produced by
# that script (run in the same directory).
# ---------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
statepoint_path = sorted(glob.glob(os.path.join(HERE, 'statepoint.*.h5')))[-1]
sp = openmc.StatePoint(statepoint_path)

material_names = {m.id: m.name for m in sp.summary.materials}


def mat_labels(material_filter):
    return [material_names[int(m)] for m in material_filter.bins]


print(f'Loaded {statepoint_path}')
print(f'k-eff: {sp.keff}\n')

## --- Reaction rate tables (absorption by material, fission in fuel) --- ##

absorption = sp.get_tally(name='absorption-by-material')
mat_filter = absorption.find_filter(openmc.MaterialFilter)
mean = absorption.get_reshaped_data(value='mean').reshape(len(mat_filter.bins), -1)[:, 0]
std = absorption.get_reshaped_data(value='std_dev').reshape(len(mat_filter.bins), -1)[:, 0]
df_abs = pd.DataFrame({
    'material': mat_labels(mat_filter),
    'absorption rate': mean,
    'std. dev.': std,
})
print('=== Absorption rate by material ===')
print(df_abs.to_string(index=False))
print()
df_abs.to_csv(os.path.join(HERE, 'absorption_by_material.csv'), index=False)

fission = sp.get_tally(name='fission-in-fuel')
mean = fission.get_reshaped_data(value='mean').reshape(-1)
std = fission.get_reshaped_data(value='std_dev').reshape(-1)
df_fis = pd.DataFrame({'score': fission.scores, 'mean': mean, 'std. dev.': std})
print('=== Fission / nu-fission rate in fuel ===')
print(df_fis.to_string(index=False))
print()
df_fis.to_csv(os.path.join(HERE, 'fission_in_fuel.csv'), index=False)

## --- Flux spectra by region --- ##

spectrum = sp.get_tally(name='flux-spectrum-by-region')
spec_mat_filter = spectrum.find_filter(openmc.MaterialFilter)
energy_filter = spectrum.find_filter(openmc.EnergyFilter)
labels = mat_labels(spec_mat_filter)
e_bins = np.array(energy_filter.bins)  # shape (n_groups, 2), eV
n_mat, n_e = len(labels), len(e_bins)

flux = spectrum.get_reshaped_data(value='mean').reshape(n_mat, n_e)

e_lo, e_hi = e_bins[:, 0], e_bins[:, 1]
e_mid = np.sqrt(np.clip(e_lo, 1e-5, None) * e_hi)  # log-midpoint
lethargy_width = np.log(e_hi / np.clip(e_lo, 1e-5, None))

fig, ax = plt.subplots(figsize=(8, 6))
for i, name in enumerate(labels):
    flux_per_lethargy = flux[i, :] / lethargy_width
    ax.step(e_mid, flux_per_lethargy, where='mid', label=name)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Flux per unit lethargy [a.u.]')
ax.set_title('Flux spectrum by region')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'flux_spectra.png'), dpi=150)
plt.close(fig)

df_spec = pd.DataFrame(flux, index=labels,
                        columns=[f'{lo:.4g}-{hi:.4g} eV' for lo, hi in e_bins])
df_spec.to_csv(os.path.join(HERE, 'flux_spectrum_by_region.csv'))

## --- Mesh flux heatmap --- ##

mesh_tally = sp.get_tally(name='flux-mesh')
mesh_filter = mesh_tally.find_filter(openmc.MeshFilter)
mesh = mesh_filter.mesh
nx, ny, nz = mesh.dimension
flux_mesh = mesh_tally.get_reshaped_data(value='mean').reshape(nx, ny, nz)

fig, ax = plt.subplots(figsize=(7, 6))
extent = (mesh.lower_left[0], mesh.upper_right[0],
          mesh.lower_left[1], mesh.upper_right[1])
# Mask cells with zero score (outside the hex boundary) so the color scale
# reflects the actual in-cell flux contrast instead of being stretched by
# the zero background.
flux_plot = np.ma.masked_equal(flux_mesh[:, :, 0].T, 0.0)
cmap = plt.get_cmap('inferno').copy()
cmap.set_bad('lightgray')
im = ax.imshow(flux_plot, origin='lower', extent=extent, cmap=cmap)
ax.set_xlabel('x [cm]')
ax.set_ylabel('y [cm]')
ax.set_title('Spatial flux distribution (mesh tally)')
fig.colorbar(im, ax=ax, label='Flux [a.u.]')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'flux_mesh_heatmap.png'), dpi=150)
plt.close(fig)

print('Wrote absorption_by_material.csv, fission_in_fuel.csv, '
      'flux_spectrum_by_region.csv, flux_spectra.png, flux_mesh_heatmap.png')
