import numpy as np
import openmc

mats = {}
mats['UN'] = openmc.Material(name='UN')
mats['UN'].add_element('U', 1.0, enrichment=19.50)
mats['UN'].add_element('N', 1.0)
temp_UN = lambda T: 14.32 / (1 + (6.9E-06+1.5E-09*(T - 700))*(T - 700))
T_UN = 900. # K
mats['UN'].temperature = T_UN
mats['UN'].set_density('g/cm3', temp_UN(mats['UN'].temperature))

print(mats['UN'])