import numpy as np
import openmc

## Materials ##
mats = {}
mats['UN'] = openmc.Material(name='UN')
mats['UN'].add_element('U', 1.0, enrichment=19.50)
mats['UN'].add_element('N', 1.0)
density_UN = lambda T: 14.32 / (1 + (6.9E-06+1.5E-09*(T - 298))*(T - 298))
T_UN = 900. # K
mats['UN'].temperature = T_UN
mats['UN'].set_density('g/cm3', density_UN(mats['UN'].temperature))

mats['YH2'] = openmc.Material(name='YH2')
mats['YH2'].add_element('Y', 1.0)
mats['YH2'].add_element('H', 2.0)
density_YH2 = lambda T: 4.085  #FIXME
T_YH2 = 900. # K
mats['YH2'].temperature = T_YH2
mats['YH2'].set_density('g/cm3', density_YH2(mats['YH2'].temperature))

mats['NaGas'] = openmc.Material(name='NaGas')
mats['NaGas'].add_element('Na', 1.0)
density_NaGas = lambda T: 0.7266  #FIXME
T_NaGas = 1200. # K
mats['NaGas'].temperature = T_NaGas
mats['NaGas'].set_density('g/cm3', density_NaGas(mats['NaGas'].temperature))

mats['ZrHx'] = openmc.Material(name='ZrHx')
mats['ZrHx'].add_element('Zr', 1.0)
mats['ZrHx'].add_element('H', 1.6)
density_ZrHx = lambda T: 5.66  #FIXME
T_ZrHx = 1000. # K
mats['ZrHx'].temperature = T_ZrHx
mats['ZrHx'].set_density('g/cm3', density_ZrHx(mats['ZrHx'].temperature))

mats['NaLiquid'] = openmc.Material(name='NaLiquid')
mats['NaLiquid'].add_element('Na', 1.0)
density_NaLiquid = lambda T: 0.970  #FIXME
T_NaLiquid = 900. # K
mats['NaLiquid'].temperature = T_NaLiquid
mats['NaLiquid'].set_density('g/cm3', density_NaLiquid(mats['NaLiquid'].temperature))

mats['Nb'] = openmc.Material(name='Nb')
mats['Nb'].add_element('Nb', 1.0)
density_Nb = lambda T: 8.57  #FIXME
T_Nb = 900. # K
mats['Nb'].temperature = T_Nb
mats['Nb'].set_density('g/cm3', density_Nb(mats['Nb'].temperature))

mats['SiC'] = openmc.Material(name='SiC')
mats['SiC'].add_element('Si', 1.0)
mats['SiC'].add_element('C' , 1.0)
density_SiC = lambda T: 3.16  #FIXME
T_SiC = 900. # K
mats['SiC'].temperature = T_SiC
mats['SiC'].set_density('g/cm3', density_SiC(mats['SiC'].temperature))

mats['W'] = openmc.Material(name='W')
mats['W'].add_element('W', 1.0)
density_W = lambda T: 19.25  #FIXME
T_W = 900. # K
mats['W'].temperature = T_W
mats['W'].set_density('g/cm3', density_W(mats['W'].temperature))

volFracs_VHTM = {
    mats['NaGas']   :(0.70**2)/1.02**2,
    mats['ZrHx']    :(0.90**2-0.70**2)/1.02**2,
    mats['NaLiquid']:(0.94**2-0.90**2)/1.02**2,
    mats['Nb']      :(0.95**2-0.94**2 + 1.01**2-1.00**2)/1.02**2,
    mats['SiC']     :(1.00**2-0.95**2)/1.02**2,
    mats['W']       :(1.02**2-1.01**2)/1.02**2,
}

materials_file = openmc.Materials(mats.values())
materials_file.cross_sections = '/home/rdwillat/openmc/XSData/endfb-viii.1-hdf5/cross_sections.xml'
# Exported after the UN fuel volume is set below, but before the geometry
# plot, which needs materials.xml on disk.

## Geometry ##
h_fuel     =   40.0
r_HPvap    = 0.70
t_wick     = 0.20
t_HPliq    = 0.04
t_liner1   = 0.01
t_envelope1= 0.05
t_liner2   = 0.02
t_fuel     = 0.75
t_envelope2= 0.40
p_pincell  = 5.5
# Radial Surfaces
s_HPvap     = openmc.ZCylinder(name='HPvap'    , r=r_HPvap)
s_HPwick    = openmc.ZCylinder(name='HPwick'   , r=r_HPvap+t_wick)
s_HPliq     = openmc.ZCylinder(name='HPliq'    , r=r_HPvap+t_wick+t_HPliq)
s_NbLiner1  = openmc.ZCylinder(name='NbLiner1' , r=r_HPvap+t_wick+t_HPliq+t_liner1)
s_envelope1 = openmc.ZCylinder(name='envelope1', r=r_HPvap+t_wick+t_HPliq+t_liner1+t_envelope1)
s_NbLiner2  = openmc.ZCylinder(name='NbLiner2' , r=r_HPvap+t_wick+t_HPliq+t_liner1+t_envelope1+t_liner1)
s_WLiner1   = openmc.ZCylinder(name='WLiner1'  , r=r_HPvap+t_wick+t_HPliq+t_liner1+t_envelope1+t_liner1+t_liner2)
s_Fuel      = openmc.ZCylinder(name='Fuel'     , r=r_HPvap+t_wick+t_HPliq+t_liner1+t_envelope1+t_liner1+t_liner2+t_fuel)
s_WLiner2   = openmc.ZCylinder(name='WLiner2'  , r=r_HPvap+t_wick+t_HPliq+t_liner1+t_envelope1+t_liner1+t_liner2+t_fuel+t_liner1)
s_NbLiner3  = openmc.ZCylinder(name='NbLiner3' , r=r_HPvap+t_wick+t_HPliq+t_liner1+t_envelope1+t_liner1+t_liner2+t_fuel+t_liner1+t_liner1)
# Axial Surfaces
s_bot = openmc.ZPlane(name='bot', z0= -h_fuel / 2.0, boundary_type='reflective')
s_top = openmc.ZPlane(name='top', z0=  h_fuel / 2.0, boundary_type='reflective')
# Planes
apothem = p_pincell / 2.0
planes_mod      = []
planes_liner    = []
planes_envelope = []
for i in range(6):
    theta = np.pi / 6.0 + i * np.pi / 3.0
    nx = np.cos(theta)
    ny = np.sin(theta)
    planes_mod+=[openmc.Plane(
        name=f'modPlane{i}',
        a=nx, b=ny, c=0.0,
        d=apothem-t_envelope2-t_liner2
    )]
    planes_liner+=[openmc.Plane(
        name=f'NbLinerPlane{i}',
        a=nx, b=ny, c=0.0,
        d=apothem-t_envelope2
    )]
    planes_envelope+=[openmc.Plane(
        name=f'envelopePlane{i}',
        a=nx, b=ny, c=0.0,
        d=apothem,
        boundary_type='reflective'
    )]
hexReg_mod = -planes_mod[0] & -planes_mod[1] & -planes_mod[2] & -planes_mod[3] & -planes_mod[4] & -planes_mod[5]
hexReg_liner = -planes_liner[0] & -planes_liner[1] & -planes_liner[2] & -planes_liner[3] & -planes_liner[4] & -planes_liner[5]
hexReg_envelope = -planes_envelope[0] & -planes_envelope[1] & -planes_envelope[2] & -planes_envelope[3] & -planes_envelope[4] & -planes_envelope[5]
# Cells
c_HPvap   = openmc.Cell(
    name  ='HPvap'    ,
    fill  =mats['NaGas']   ,
    region= -s_HPvap
          & +s_bot & -s_top
)
c_HPwick  = openmc.Cell(
    name  ='HPwick'   ,
    fill  =mats['ZrHx']   ,
    region= -s_HPwick & +s_HPvap
          & +s_bot & -s_top
)
c_HPliq   = openmc.Cell(
    name  ='HPliq'    ,
    fill  =mats['NaLiquid']   ,
    region= -s_HPliq & +s_HPwick
          & +s_bot & -s_top
)
c_NbLiner1 = openmc.Cell(
    name  ='NbLiner1' ,
    fill  =mats['Nb'] ,
    region=-s_NbLiner1 & +s_HPliq
          & +s_bot & -s_top
)
c_SiC_envelope1 = openmc.Cell(
    name  ='SiC_envelope1' ,
    fill  =mats['SiC'] ,
    region= -s_envelope1 & +s_NbLiner1
          & +s_bot & -s_top
)
c_NbLiner2 = openmc.Cell(
    name  ='NbLiner2' ,
    fill  =mats['Nb'] ,
    region= -s_NbLiner2 & +s_envelope1
          & +s_bot & -s_top
)
c_WLiner1 = openmc.Cell(
    name  ='WLiner1' ,
    fill  =mats['W'] ,
    region= -s_WLiner1 & +s_NbLiner2
          & +s_bot & -s_top
)
c_Fuel = openmc.Cell(
    name  ='Fuel' ,
    fill  =mats['UN'] ,
    region= -s_Fuel & +s_WLiner1
          & +s_bot & -s_top
)
# Exact UN fuel volume (annulus between s_WLiner1 and s_Fuel, full h_fuel
# height), needed for the specific-power postprocessing metric.
mats['UN'].volume = np.pi * (s_Fuel.r ** 2 - s_WLiner1.r ** 2) * h_fuel
c_WLiner2 = openmc.Cell(
    name  ='WLiner2' ,
    fill  =mats['W'] ,
    region= -s_WLiner2 & +s_Fuel
          & +s_bot & -s_top
)
c_NbLiner3 = openmc.Cell(
    name  ='NbLiner3' ,
    fill  =mats['Nb'] ,
    region= -s_NbLiner3 & +s_WLiner2
          & +s_bot & -s_top
)
c_mod = openmc.Cell(
    name  ='mod' ,
    fill  =mats['YH2'] ,
    region= hexReg_mod & +s_NbLiner3
          & +s_bot & -s_top
)
c_NbLiner4 = openmc.Cell(
    name  ='NbLiner4' ,
    fill  =mats['Nb'] ,
    region= hexReg_liner & ~hexReg_mod
          & +s_bot & -s_top
)
c_SiC_envelope2 = openmc.Cell(
    name  ='SiC_envelope2' ,
    fill  =mats['SiC'] ,
    region= hexReg_envelope & ~hexReg_liner
          & +s_bot & -s_top
)
cells = [
    c_HPvap,
    c_HPwick,
    c_HPliq,
    c_NbLiner1,
    c_SiC_envelope1,
    c_NbLiner2,
    c_WLiner1,
    c_Fuel,
    c_WLiner2,
    c_NbLiner3,
    c_mod,
    c_NbLiner4,
    c_SiC_envelope2,
]
# Universes
pin_cell_universe = openmc.Universe(
    name='pin_cell',
    universe_id=0,
    cells=cells,
)
geometry_file = openmc.Geometry(root=pin_cell_universe)
geometry_file.export_to_xml()
materials_file.export_to_xml()

## Plot ##
plot = openmc.Plot()
plot.filename = 'pincell'
plot.width = (p_pincell * 1.1, p_pincell * 1.1)
plot.pixels = (1200, 1200)
plot.color_by = 'material'
plot.to_ipython_image()
openmc.Plots([plot]).export_to_xml()

## Settings ##
settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.particles = 5000
settings.batches = 110
settings.inactive = 10
settings.temperature = {'default': 900.0, 'method': 'interpolation'}
settings.export_to_xml()

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

# Same 11-group energy structure (MeV -> eV converted) as the pinCell cases,
# for direct comparison
energy_bins = [1.0e-5, 8.00e-2, 1.80e-1, 6.25e-1, 1.30e+0, 4.00e+0,
               1.4873e+2, 9.118e+3, 1.83e+5, 5.00e+5, 1.353e+6, 2.0e+7]
energy_filter = openmc.EnergyFilter(energy_bins)
spectrum_materials = [mats['UN'], mats['YH2'], mats['NaGas'], mats['NaLiquid'],
                      mats['ZrHx']]
mat_filter_spectrum = openmc.MaterialFilter(spectrum_materials)
tally_spectrum = openmc.Tally(name='flux-spectrum-by-region')
tally_spectrum.filters = [mat_filter_spectrum, energy_filter]
tally_spectrum.scores = ['flux']
tallies.append(tally_spectrum)

mesh = openmc.RegularMesh()
mesh.dimension = (200, 200, 1)
mesh.lower_left = (-p_pincell / 2 * 1.05, -p_pincell / 2 * 1.05, -h_fuel / 2)
mesh.upper_right = (p_pincell / 2 * 1.05, p_pincell / 2 * 1.05, h_fuel / 2)
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
