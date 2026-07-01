import numpy as np
import openmc

# Approximate pin-cell replica of SerpentTests/HPMR.sss: a triangular-pitch
# "flower" unit cell (1 heat pipe pin + ring of 6 alternating fuel/moderator
# pins, pitch 2.3 cm) with reflective boundaries, for a k-inf check -- not
# the full core. TRISO particles are stochastically re-packed at PF=40%
# rather than reusing Serpent's exact coordinate file.

## Materials ##
mats = {}

mats['fuel'] = openmc.Material(name='fuel')
mats['fuel'].set_density('g/cm3', 10.7440)
mats['fuel'].add_nuclide('U235', 6.8794e-02, 'ao')
mats['fuel'].add_nuclide('U238', 2.7604e-01, 'ao')
mats['fuel'].add_nuclide('C12', 1.3793e-01, 'ao')
mats['fuel'].add_nuclide('O16', 5.1724e-01, 'ao')
mats['fuel'].temperature = 900.0

mats['buffer'] = openmc.Material(name='buffer')
mats['buffer'].set_density('g/cm3', 1.0400)
mats['buffer'].add_nuclide('C12', 1.0, 'ao')
mats['buffer'].temperature = 900.0

mats['PyC1'] = openmc.Material(name='PyC1')
mats['PyC1'].set_density('g/cm3', 1.8820)
mats['PyC1'].add_nuclide('C12', 1.0, 'ao')
mats['PyC1'].temperature = 900.0

mats['SiC'] = openmc.Material(name='SiC')
mats['SiC'].set_density('g/cm3', 3.1710)
mats['SiC'].add_nuclide('Si28', 0.4611, 'ao')
mats['SiC'].add_nuclide('Si29', 0.0234, 'ao')
mats['SiC'].add_nuclide('Si30', 0.0154, 'ao')
mats['SiC'].add_nuclide('C12', 0.5, 'ao')
mats['SiC'].temperature = 900.0

mats['PyC2'] = openmc.Material(name='PyC2')
mats['PyC2'].set_density('g/cm3', 1.8820)
mats['PyC2'].add_nuclide('C12', 1.0, 'ao')
mats['PyC2'].temperature = 900.0

# Graphite background inside the TRISO compact
mats['matrix_pin'] = openmc.Material(name='matrix_pin')
mats['matrix_pin'].set_density('g/cm3', 1.8060)
mats['matrix_pin'].add_nuclide('C12', 0.9999997, 'ao')
mats['matrix_pin'].add_nuclide('B10', 0.0000003, 'ao')
mats['matrix_pin'].add_s_alpha_beta('c_Graphite')
mats['matrix_pin'].temperature = 900.0

# Graphite filling the rest of each pin cell outside the compact
mats['matrix'] = openmc.Material(name='matrix')
mats['matrix'].set_density('g/cm3', 1.8060)
mats['matrix'].add_nuclide('C12', 0.9999997, 'ao')
mats['matrix'].add_nuclide('B10', 0.0000003, 'ao')
mats['matrix'].add_s_alpha_beta('c_Graphite')
mats['matrix'].temperature = 700.0

mats['moderator'] = openmc.Material(name='moderator')
mats['moderator'].set_density('g/cm3', 4.0850)
mats['moderator'].add_nuclide('Y89', 0.357142857, 'ao')
mats['moderator'].add_nuclide('H1', 0.642857143, 'ao')
mats['moderator'].add_s_alpha_beta('c_H_in_YH2')
mats['moderator'].add_s_alpha_beta('c_Y_in_YH2')
mats['moderator'].temperature = 700.0

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

# Smeared He gas + stainless steel annulus around the moderator pin
mats['shell_air_mod'] = _he_steel_mix('shell_air_mod', 2.290e-02, [
    ('He4', 1.983e-05), ('C12', 4.349e-05),
    ('Si28', 2.121e-04), ('Si29', 1.081e-05), ('Si30', 7.130e-06),
    ('P31', 9.454e-06),
    ('S32', 5.653e-06), ('S33', 4.464e-08), ('S34', 2.530e-07), ('S36', 6.867e-10),
    ('Cr50', 1.810e-04), ('Cr52', 3.491e-03), ('Cr53', 3.958e-04), ('Cr54', 9.852e-05),
    ('Mn55', 2.352e-04),
    ('Fe54', 8.929e-04), ('Fe56', 1.400e-02), ('Fe57', 3.236e-04), ('Fe58', 3.053e-05),
    ('Ni58', 1.773e-03), ('Ni60', 6.831e-04), ('Ni61', 2.970e-05), ('Ni62', 9.469e-05), ('Ni64', 2.412e-05),
    ('Mo92', 4.864e-05), ('Mo94', 3.043e-05), ('Mo95', 5.269e-05), ('Mo96', 5.535e-05),
    ('Mo97', 3.181e-05), ('Mo98', 8.065e-05), ('Mo100', 3.234e-05),
])

materials = openmc.Materials(mats.values())
materials.cross_sections = '/home/rdwillat/openmc/XSData/endfb-viii.1-hdf5/cross_sections.xml'
# Exported after the fuel volume is set below, but before the geometry plot,
# which needs materials.xml on disk.

## Geometry ##

pitch = 2.3     # cm, triangular pin pitch
h_pin = 4.0     # cm, axial slice height of the reflective unit cell
pf = 0.40       # TRISO packing fraction

# TRISO layer radii, cm
r_kernel = 0.02125
r_buffer = 0.03125
r_PyC1 = 0.03525
r_SiC = 0.03875
r_PyC2 = 0.04275

r_compact = 1.00
r_mod_core = 0.825
r_mod_shell = 0.92
r_hp_core = 0.97
r_hp_shell = 1.07

z_lo = openmc.ZPlane(z0=-h_pin / 2, boundary_type='reflective')
z_hi = openmc.ZPlane(z0=h_pin / 2, boundary_type='reflective')

# TRISO particle: concentric fuel/buffer/PyC1/SiC/PyC2 spheres
s_kernel = openmc.Sphere(r=r_kernel)
s_buffer = openmc.Sphere(r=r_buffer)
s_PyC1 = openmc.Sphere(r=r_PyC1)
s_SiC = openmc.Sphere(r=r_SiC)
s_PyC2 = openmc.Sphere(r=r_PyC2)

triso_universe = openmc.Universe(cells=[
    openmc.Cell(fill=mats['fuel'], region=-s_kernel),
    openmc.Cell(fill=mats['buffer'], region=+s_kernel & -s_buffer),
    openmc.Cell(fill=mats['PyC1'], region=+s_buffer & -s_PyC1),
    openmc.Cell(fill=mats['SiC'], region=+s_PyC1 & -s_SiC),
    openmc.Cell(fill=mats['PyC2'], region=+s_SiC & -s_PyC2),
])

s_compact = openmc.ZCylinder(r=r_compact)

# Pack spheres shrunk by r_PyC2 on every side so no TRISO particle's outer
# surface extends past the compact (otherwise OpenMC warns of a particle
# outside the lattice).
s_pack = openmc.ZCylinder(r=r_compact - r_PyC2)
z_pack_lo = openmc.ZPlane(z0=-h_pin / 2 + r_PyC2)
z_pack_hi = openmc.ZPlane(z0=h_pin / 2 - r_PyC2)
pack_region = -s_pack & +z_pack_lo & -z_pack_hi
triso_centers = openmc.model.pack_spheres(radius=r_PyC2, region=pack_region, pf=pf)
trisos = [openmc.model.TRISO(r_PyC2, triso_universe, c) for c in triso_centers]

# fuel_pin is instanced 3 times in the hex flower below, and material-filter
# tallies sum over all instances, so the fuel volume must count every pin.
n_fuel_pins = 3  # positions 0, 2, 4 of the 6-pin outer ring (see hex_lattice)
n_kernels = len(triso_centers)
mats['fuel'].volume = n_fuel_pins * n_kernels * (4.0 / 3.0) * np.pi * r_kernel ** 3

n_xy = 10
n_z = max(1, round(h_pin / (2 * r_compact / n_xy)))
triso_lattice = openmc.model.create_triso_lattice(
    trisos,
    lower_left=(-r_compact, -r_compact, -h_pin / 2),
    pitch=(2 * r_compact / n_xy, 2 * r_compact / n_xy, h_pin / n_z),
    shape=(n_xy, n_xy, n_z),
    background=mats['matrix_pin'],
)
triso_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=mats['matrix_pin'])])

# Fuel compact pin
fuel_pin = openmc.Universe(cells=[
    openmc.Cell(region=-s_compact, fill=triso_lattice),
    openmc.Cell(region=+s_compact, fill=mats['matrix']),
])

# Moderator pin
s_mod_core = openmc.ZCylinder(r=r_mod_core)
s_mod_shell = openmc.ZCylinder(r=r_mod_shell)
mod_pin = openmc.Universe(cells=[
    openmc.Cell(region=-s_mod_core, fill=mats['moderator']),
    openmc.Cell(region=+s_mod_core & -s_mod_shell, fill=mats['shell_air_mod']),
    openmc.Cell(region=+s_mod_shell, fill=mats['matrix']),
])

# Heat pipe pin
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

mat_filter_fuel = openmc.MaterialFilter([mats['fuel']])
tally_fission = openmc.Tally(name='fission-in-fuel')
tally_fission.filters = [mat_filter_fuel]
tally_fission.scores = ['fission', 'nu-fission', 'kappa-fission']
tallies.append(tally_fission)

tally_heating = openmc.Tally(name='heating-by-material')
tally_heating.filters = [mat_filter_all]
tally_heating.scores = ['heating-local']
tallies.append(tally_heating)

# 11-group structure from Serpent's `set nfg` card (MeV -> eV, with 1e-5 eV /
# 20 MeV as the overall bounds)
energy_bins = [1.0e-5, 8.00e-2, 1.80e-1, 6.25e-1, 1.30e+0, 4.00e+0,
               1.4873e+2, 9.118e+3, 1.83e+5, 5.00e+5, 1.353e+6, 2.0e+7]
energy_filter = openmc.EnergyFilter(energy_bins)
spectrum_materials = [mats['fuel'], mats['moderator'], mats['hp_vp_liq_wick'],
                      mats['matrix'], mats['matrix_pin']]
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

# openmc.run()
