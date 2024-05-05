import openmc
import math

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

openmc.Materials.cross_sections='/home/shared/endfb71_hdf5/cross_sections.xml'
# Create the pin cell geometry


# Create the pin_cell materials
fuel = openmc.Material(name='UO2')
fuel.set_density('g/cm3', 10.29769)
fuel.add_element('U', 1.0, enrichment=4.5)
fuel.add_element('O', 2.0)

clad = openmc.Material(name='Zircaloy')
clad.set_density('g/cm3', 6.55)
clad.add_nuclide('Zr90', 2.1827e-2)
clad.add_nuclide('Zr91', 4.7600e-3)
clad.add_nuclide('Zr92', 7.2758e-3)
clad.add_nuclide('Zr94', 7.3734e-3)
clad.add_nuclide('Zr96', 1.1879e-3)

hot_water = openmc.Material(name='Hot borated water')
hot_water.set_density('g/cm3', 0.740582)
hot_water.add_nuclide('H1', 4.9457e-2)
hot_water.add_nuclide('O16', 2.4672e-2)
hot_water.add_nuclide('B10', 8.0042e-6)
hot_water.add_nuclide('B11', 3.2218e-5)
hot_water.add_s_alpha_beta('c_H_in_H2O')



# Instantiate ZCylinder surfaces
pitch = 1.26
fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.40947, name='Fuel OR')
clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR')
left = openmc.XPlane(x0=-pitch/2, name='left', boundary_type='reflective')
right = openmc.XPlane(x0=pitch/2, name='right', boundary_type='reflective')
bottom = openmc.YPlane(y0=-pitch/2, name='bottom',
                       boundary_type='reflective')
top = openmc.YPlane(y0=pitch/2, name='top', boundary_type='reflective')



# Instantiate Cells
fuel_pin = openmc.Cell(name='Fuel', fill=fuel)
cladding = openmc.Cell(name='Cladding', fill=clad)
water = openmc.Cell(name='Water', fill=hot_water)



# Use surface half-spaces to define regions
fuel_pin.region = -fuel_or
cladding.region = +fuel_or & -clad_or
water.region = +clad_or & +left & -right & +bottom & -top


# Create the model

pin_cell=openmc.model.Model()

# Define the materials file
pin_cell.materials=(fuel, clad, hot_water)

# Create root universe
pin_cell.geometry.root_universe = openmc.Universe(0, name='root universe')
pin_cell.geometry.root_universe.add_cells([fuel_pin, cladding, water])



# Settings
pin_cell.settings.batches = 200
pin_cell.settings.inactive = 20
pin_cell.settings.particles = 5000
pin_cell.settings.source = openmc.Source(space=openmc.stats.Box(
    [-pitch/2, -pitch/2, -1], [pitch/2, pitch/2, 1], only_fissionable=True))

plot = openmc.Plot.from_geometry(pin_cell.geometry)
plot.pixels = (300, 300)
plot.color_by = 'material'
pin_cell.plots.append(plot)


# Export to xml

pin_cell.materials.export_to_xml()
pin_cell.geometry.export_to_xml()
pin_cell.settings.export_to_xml()

# Plot
pin_cell.geometry.root_universe.plot(width=(-1.0, 1.0))

# Set the volume of the depletable material 

rod_radius=0.40947 # Value from AP1000 documentation
pin_cell.materials[0].volume=(math.pi*rod_radius**2)*pitch # 3D representation with the pitch of the fuel pin
volume=pin_cell.materials[0].volume
print("Volume HM :", volume)


##############################################################################################################
                                          # DEPLETION 
##############################################################################################################

import openmc.deplete

chain_file="/home/chloec/Sources/chain_endfb71_pwr.xml"

openmc.deplete.Chain.from_xml(chain_file)

operator=openmc.deplete.Operator(pin_cell, chain_file)

# Set the power schedule for depletion
# values where taken from AP1000 documentation

# The linear power of AP1000 is 5.72 kW/ft ie 187.6 W/cm

linear_power=5.72*1000/30.48 #W/cm 
power=linear_power*pitch

BU_final=50 #MWd/kgHM
N_days=(BU_final*operator.heavy_metal*1000)/power

end_step=N_days-1219


time_operation = ([3/24]*4)+([4/24]*3)+([6/24]*4)+([12/24]*2)+([1])+([3]*1)+([7]*3)+([14]*6)+([28]*8)+([29])+\
             ([31])+([31])+ ([31])+([30])+([31])+([30])+([31])+([31])+([30])+([31]) +\
             ([59])+([61])+([61])+([62])+([61])+([61])+\
             ([60])+([61])+ ([61])+ ([end_step])

time_cooling=[0.04, 1.0, 7.0, 30.0, 30.0, 300.0]+([365]*4)+\
            [36.5, 146.0, 91.25, 91.25, 730, 730, 547.5, 1277.5, 1825, 1095, 730, 1825, 1825, 1825, 1825, 1825]+\
            [3650, 3650, 3650, 3650, 3650]+\
            [9125, 9125, 9125, 9125, 9125, 9125]+\
            [18250, 18250, 18250, 18250, 18250]+\
            [36500, 36500, 36500, 36500, 36500]+\
            [91250, 91250, 91250, 91250, 91250, 91250] +\
            [182500, 365000, 365000, 365000, 365000, 365000, 365000, 365000, 3650000]+\
            [10950000, 18250000, 36500000, 109500000, 182500000]
time_steps=time_operation+time_cooling


power_operation=[power]*len(time_operation)
power_discharge=[0]*len(time_cooling)
power_reactor=power_operation+power_discharge



integrator=openmc.deplete.CECMIntegrator(operator, time_steps, power_reactor, timestep_units='d')
print("Linear power: " + str(linear_power))
print("Running time: "+str(N_days) + ' nombre de jours simulations '+str (sum(time_operation)))
print("operator.heavy_metal : " + str(operator.heavy_metal))

openmc.run()
#integrator.integrate()



