import numpy as np
import openmc

# Variant of HPMR_pinCell/HPMR_pincell.py: same hex "flower" unit cell, but
# with the TRISO compact replaced by a solid UN fuel pellet (composition from
# HPMR_INFINIT.py) clad in SiC, and the fuel/moderator pins' He+steel smeared
# shells replaced by SiC cladding. Heat pipe pin, matrix, lattice, and tally
# set are unchanged so the two cases are directly comparable.

## Materials ##
mats = {}

# Solid UN fuel pellet (composition/density as in HPMR_INFINIT.py)
mats['UN'] = openmc.Material(name='UN')
mats['UN'].add_element('U', 1.0, enrichment=19.50)
mats['UN'].add_element('N', 1.0)
density_UN = lambda T: 14.32 / (1 + (6.9E-06 + 1.5E-09 * (T - 298)) * (T - 298))
mats['UN'].temperature = 900.0
mats['UN'].set_density('g/cm3', density_UN(mats['UN'].temperature))

# SiC cladding (same composition as the TRISO SiC layer in HPMR_pincell.py).
# Separate instances per pin type so by-material tallies still distinguish
# fuel-clad from moderator-clad.
def _sic(name, temperature):
    m = openmc.Material(name=name)
    m.set_density('g/cm3', 3.1710)
    m.add_nuclide('Si28', 0.4611, 'ao')
    m.add_nuclide('Si29', 0.0234, 'ao')
    m.add_nuclide('Si30', 0.0154, 'ao')
    m.add_nuclide('C12', 0.5, 'ao')
    m.temperature = temperature
    return m

mats['SiC_fuel_clad'] = _sic('SiC_fuel_clad', 900.0)
mats['SiC_mod_clad'] = _sic('SiC_mod_clad', 700.0)

mats['moderator'] = openmc.Material(name='moderator')
mats['moderator'].set_density('g/cm3', 4.0850)
mats['moderator'].add_nuclide('Y89', 0.357142857, 'ao')
mats['moderator'].add_nuclide('H1', 0.642857143, 'ao')
mats['moderator'].add_s_alpha_beta('c_H_in_YH2')
mats['moderator'].add_s_alpha_beta('c_Y_in_YH2')
mats['moderator'].temperature = 700.0

mats['matrix'] = openmc.Material(name='matrix')
mats['matrix'].set_density('g/cm3', 1.8060)
mats['matrix'].add_nuclide('C12', 0.9999997, 'ao')
mats['matrix'].add_nuclide('B10', 0.0000003, 'ao')
mats['matrix'].add_s_alpha_beta('c_Graphite')
mats['matrix'].temperature = 700.0

# Heat pipe working fluid (Na vapor/liquid + wick), homogenized
mats['hp_vp_liq_wick'] = openmc.Material(name='hp_vp_liq_wick')
mats['hp_vp_liq_wick'].set_density('atom/b-cm', 8.324e-03)
for nuclide, frac in [
    ('C12', 3.808e-06), ('Si28', 1.856e-05), ('Si29', 9.463e-07), ('Si30', 6.242e-07),
    ('P31', 8.277e-07), ('S32', 4.949e-07), ('S33', 3.908e-09), ('S34', 2.214e-08),
    ('S36', 6.012e-11), ('Cr50', 1.585e-05), ('Cr52', 3.055e-04), ('Cr53', 3.465e-05),
    ('Cr54', 8.625e-06), ('Mn55', 2.059e-05), ('Fe54', 7.814e-05), ('Fe56', 1.226e-03),
    ('Fe57', 2.833e-05), ('Fe58', 2.673e-06), ('Ni58', 1.553e-04), ('Ni60', 5.980e-05),
    ('Ni61', 2.600e-06), ('Ni62', 8.291e-06), ('Ni64', 2.112e-06), ('Mo92', 4.258e-06),
    ('Mo94', 2.664e-06), ('Mo95', 4.613e-06), ('Mo96', 4.845e-06), ('Mo97', 2.784e-06),
    ('Mo98', 7.058e-06), ('Mo100', 2.831e-06), ('K39', 5.895e-03), ('K40', 7.586e-07),
    ('K41', 4.254e-04),
]:
    mats['hp_vp_liq_wick'].add_nuclide(nuclide, frac, 'ao')
mats['hp_vp_liq_wick'].temperature = 700.0


def _he_steel_mix(name, density, composition):
    m = openmc.Material(name=name)
    m.set_density('atom/b-cm', density)
    for nuclide, frac in composition:
        m.add_nuclide(nuclide, frac, 'ao')
    m.temperature = 700.0
    return m


# Smeared He gas + stainless steel annulus around the heat pipe
mats['shell_air_hp'] = _he_steel_mix('shell_air_hp', 6.771e-02, [
    ('He4', 5.629e-06), ('C12', 1.287e-04),
    ('Si28', 6.276e-04), ('Si29', 3.199e-05), ('Si30', 2.110e-05),
    ('P31', 2.798e-05),
    ('S32', 1.673e-05), ('S33', 1.321e-07), ('S34', 7.486e-07), ('S36', 2.032e-09),
    ('Cr50', 5.357e-04), ('Cr52', 1.033e-02), ('Cr53', 1.171e-03), ('Cr54', 2.916e-04),
    ('Mn55', 6.960e-04),
    ('Fe54', 2.643e-03), ('Fe56', 4.145e-02), ('Fe57', 9.576e-04), ('Fe58', 9.034e-05),
    ('Ni58', 5.248e-03), ('Ni60', 2.022e-03), ('Ni61', 8.789e-05), ('Ni62', 2.802e-04), ('Ni64', 7.139e-05),
    ('Mo92', 1.439e-04), ('Mo94', 9.004e-05), ('Mo95', 1.559e-04), ('Mo96', 1.638e-04),
    ('Mo97', 9.413e-05), ('Mo98', 2.387e-04), ('Mo100', 9.570e-05),
])

materials = openmc.Materials(mats.values())
materials.cross_sections = '/home/rdwillat/openmc/XSData/endfb-viii.1-hdf5/cross_sections.xml'
# Exported after the UN fuel volume is set below, but before the geometry
# plot, which needs materials.xml on disk.

## Geometry ##

pitch = 2.3     # cm, triangular pin pitch
h_pin = 4.0     # cm, axial slice height of the reflective unit cell

# Fuel pin: solid UN pellet clad in SiC. Clad outer radius matches the old
# TRISO compact envelope (r_compact = 1.00 cm) so the pin keeps the same
# footprint as HPMR_pincell.py's fuel_pin.
r_fuel_clad = 1.00
t_fuel_clad = 0.05
r_fuel_core = r_fuel_clad - t_fuel_clad

# Moderator pin: same core/shell radii as HPMR_pincell.py, shell swapped from
# He+steel smear to SiC.
r_mod_core = 0.825
r_mod_clad = 0.92

# Heat pipe pin: unchanged from HPMR_pincell.py.
r_hp_core = 0.97
r_hp_shell = 1.07

z_lo = openmc.ZPlane(z0=-h_pin / 2, boundary_type='reflective')
z_hi = openmc.ZPlane(z0=h_pin / 2, boundary_type='reflective')

s_fuel_core = openmc.ZCylinder(r=r_fuel_core)
s_fuel_clad = openmc.ZCylinder(r=r_fuel_clad)
fuel_pin = openmc.Universe(cells=[
    openmc.Cell(region=-s_fuel_core, fill=mats['UN']),
    openmc.Cell(region=+s_fuel_core & -s_fuel_clad, fill=mats['SiC_fuel_clad']),
    openmc.Cell(region=+s_fuel_clad, fill=mats['matrix']),
])

# fuel_pin is instanced 3 times in the hex flower below, and material-filter
# tallies sum over all instances, so the fuel volume must count every pin.
n_fuel_pins = 3  # positions 0, 2, 4 of the 6-pin outer ring (see hex_lattice)
mats['UN'].volume = n_fuel_pins * np.pi * r_fuel_core ** 2 * h_pin

s_mod_core = openmc.ZCylinder(r=r_mod_core)
s_mod_clad = openmc.ZCylinder(r=r_mod_clad)
mod_pin = openmc.Universe(cells=[
    openmc.Cell(region=-s_mod_core, fill=mats['moderator']),
    openmc.Cell(region=+s_mod_core & -s_mod_clad, fill=mats['SiC_mod_clad']),
    openmc.Cell(region=+s_mod_clad, fill=mats['matrix']),
])

s_hp_core = openmc.ZCylinder(r=r_hp_core)
s_hp_shell = openmc.ZCylinder(r=r_hp_shell)
hp_pin = openmc.Universe(cells=[
    openmc.Cell(region=-s_hp_core, fill=mats['hp_vp_liq_wick']),
    openmc.Cell(region=+s_hp_core & -s_hp_shell, fill=mats['shell_air_hp']),
    openmc.Cell(region=+s_hp_shell, fill=mats['matrix']),
])

# Hex "flower" unit cell: 1 heat pipe center, ring of 6 alternating
# fuel/moderator pins, triangular pitch
num_rings = 2  # center + 1 ring of neighbors
hex_lattice = openmc.HexLattice()
hex_lattice.orientation = 'y'
hex_lattice.center = (0., 0.)
hex_lattice.pitch = (pitch,)
hex_lattice.universes = [
    [fuel_pin, mod_pin, fuel_pin, mod_pin, fuel_pin, mod_pin],
    [hp_pin],
]
hex_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=mats['matrix'])])

# Flat-to-flat width of an N-ring hex cluster of the given pitch is
# (2*N - 1) * pitch; edge_length = flat_to_flat / sqrt(3) for orientation='y'.
flat_to_flat = (2 * num_rings - 1) * pitch
edge_length = flat_to_flat / np.sqrt(3)
outer_boundary = -openmc.model.HexagonalPrism(
    edge_length=edge_length, orientation='y', boundary_type='reflective')

root_cell = openmc.Cell(fill=hex_lattice, region=outer_boundary & +z_lo & -z_hi)
geometry = openmc.Geometry(openmc.Universe(cells=[root_cell]))
geometry.export_to_xml()
materials.export_to_xml()

## Settings ##
settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.particles = 5000
settings.batches = 110
settings.inactive = 10
settings.temperature = {'default': 700.0, 'method': 'interpolation'}
settings.export_to_xml()

## Plot ##
plot = openmc.Plot()
plot.filename = 'pincell'
plot.width = (flat_to_flat * 1.1, flat_to_flat * 1.1)
plot.pixels = (1200, 1200)
plot.color_by = 'material'
plot.to_ipython_image()
openmc.Plots([plot]).export_to_xml()

## Tallies ##
tallies = openmc.Tallies()

mat_filter_all = openmc.MaterialFilter(list(mats.values()))
tally_absorption = openmc.Tally(name='absorption-by-material')
tally_absorption.filters = [mat_filter_all]
tally_absorption.scores = ['absorption']
tallies.append(tally_absorption)

mat_filter_fuel = openmc.MaterialFilter([mats['UN']])
tally_fission = openmc.Tally(name='fission-in-fuel')
tally_fission.filters = [mat_filter_fuel]
tally_fission.scores = ['fission', 'nu-fission', 'kappa-fission']
tallies.append(tally_fission)

tally_heating = openmc.Tally(name='heating-by-material')
tally_heating.filters = [mat_filter_all]
tally_heating.scores = ['heating-local']
tallies.append(tally_heating)

# Same 11-group energy structure as HPMR_pincell.py/HPMR_INFINIT.py, for
# direct comparison
energy_bins = [1.0e-5, 8.00e-2, 1.80e-1, 6.25e-1, 1.30e+0, 4.00e+0,
               1.4873e+2, 9.118e+3, 1.83e+5, 5.00e+5, 1.353e+6, 2.0e+7]
energy_filter = openmc.EnergyFilter(energy_bins)
spectrum_materials = [mats['UN'], mats['moderator'], mats['hp_vp_liq_wick'],
                      mats['matrix'], mats['SiC_fuel_clad']]
mat_filter_spectrum = openmc.MaterialFilter(spectrum_materials)
tally_spectrum = openmc.Tally(name='flux-spectrum-by-region')
tally_spectrum.filters = [mat_filter_spectrum, energy_filter]
tally_spectrum.scores = ['flux']
tallies.append(tally_spectrum)

mesh = openmc.RegularMesh()
mesh.dimension = (200, 200, 1)
mesh.lower_left = (-flat_to_flat / 2 * 1.05, -flat_to_flat / 2 * 1.05, -h_pin / 2)
mesh.upper_right = (flat_to_flat / 2 * 1.05, flat_to_flat / 2 * 1.05, h_pin / 2)
mesh_filter = openmc.MeshFilter(mesh)
tally_mesh_flux = openmc.Tally(name='flux-mesh')
tally_mesh_flux.filters = [mesh_filter]
tally_mesh_flux.scores = ['flux']
tallies.append(tally_mesh_flux)

tally_mesh_heating = openmc.Tally(name='heating-mesh')
tally_mesh_heating.filters = [mesh_filter]
tally_mesh_heating.scores = ['heating-local']
tallies.append(tally_mesh_heating)

tallies.export_to_xml()
