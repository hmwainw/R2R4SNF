import openmc
import numpy as np
import time

chain_file = "/Users/krishnasunder/Desktop/workspace/R2R4SNF/nuclear_data/chain_endfb71_sfr.xml"
#openmc.Materials.cross_sections=

# Temperature definitions

## Temperature of the active core
fueltemp = 582.5+273.15 #K
structemp=735.65 #K
cooltemp = 705.65 #K

refltemp = 705.65 #K
abstemp = 705.65 #K
shieldtemp= 628.15 #K

R = 8.314 #m3Pa /K/mol
T = cooltemp
M = 22.98 #g/mol
P = 6.39e6 #Pa
rho = M * P / (R*T) / 1e6 #g/cm3
print(rho)


# Materials definitions

fuel=openmc.Material(name='UPuZr') # Select zone 20% Pu enriched in the ABTR benchmark model
fuel.add_nuclide('U234', 0.0000005)
fuel.add_nuclide('U235', 0.0009796)
fuel.add_nuclide('U236', 0.0000573)
fuel.add_nuclide('U238', 0.5957002)
fuel.add_nuclide('Np237', 0.0000964)
fuel.add_nuclide('Pu238', 0.0000305)
fuel.add_nuclide('Pu239', 0.1328137)
fuel.add_nuclide('Pu240', 0.0140674)
fuel.add_nuclide('Pu241', 0.0009111)
fuel.add_nuclide('Pu242', 0.0000636)
fuel.add_nuclide('Am241', 0.0000687)
fuel.add_nuclide('Am242', 0.0000013)
fuel.add_nuclide('Am243', 0.0000023)
fuel.add_nuclide('Cm242', 0.0000017)
fuel.add_nuclide('Cm244', 0.0000002)
fuel.add_element('Zr', 0.2249272)
fuel.add_element('Mo', 0.0302783)
fuel.set_density('g/cm3', 15.8)


HT_9=openmc.Material(name='HT_9')
HT_9.add_element('Fe', 0.8453, 'wo')
HT_9.add_element('Cr', 0.122, 'wo')
HT_9.add_element('Mo', 0.0090, 'wo')
HT_9.add_element('W', 0.0050, 'wo')
HT_9.add_element('Ni', 0.0069, 'wo')
HT_9.add_element('V', 0.0029, 'wo')
HT_9.add_element('C', 0.0016, 'wo')
HT_9.add_element('Si', 0.0004, 'wo')
HT_9.add_element('Mn', 0.0058, 'wo')
HT_9.add_element('N', 0.0011, 'wo')
HT_9.set_density('g/cm3', 6.551)
HT_9.temperature=structemp

reflector=openmc.Material(name='reflector')
reflector.add_element('Fe', 0.85582)
reflector.add_element('Ni', 0.00528)
reflector.add_element('Cr', 0.12725)
reflector.add_nuclide('Mn55',0.00564)
reflector.add_element('Mo', 0.00602)
reflector.set_density('g/cm3', 6.551)
reflector.temperature=refltemp

absorber=openmc.Material(name='absorber')
absorber.add_element('C', 0.0469)
absorber.add_nuclide('B10', 0.190)
absorber.add_nuclide('B11', 0.763)
absorber.temperature=abstemp


ref_homo=openmc.Material(name='reflector_homo')
ref_homo.add_element('Na', 0.048381296)
ref_homo.add_element('Fe', 0.814414793)
ref_homo.add_element('Ni', 0.005021413)
ref_homo.add_element('Cr', 0.12109295)
ref_homo.add_nuclide('Mn55',0.005364603)
ref_homo.add_element('Mo', 0.005724945)
ref_homo.temperature=refltemp

coolant=openmc.Material(name='coolant')
coolant.add_element('Na', 1.0)
coolant.set_density('g/cm3', rho )
coolant.temperature=cooltemp

shield=openmc.Material(name='shield')
shield.add_element('C', 0.02706)
shield.add_nuclide('B10', 0.1141)
shield.add_nuclide('B11', 0.4593)
shield.add_element('Na', 0.0537)
shield.add_element('Fe', 0.2959)
shield.add_element('Ni', 0.0018)
shield.add_element('Cr', 0.0440)
shield.add_element('Mo', 0.0020)
shield.add_nuclide('Mn55', 0.0019)
shield.set_density('g/cm3', 4.14)
shield.temperature=shieldtemp

barrel=openmc.Material(name='barrel')
barrel.add_element('Na',0.5617 )
barrel.add_element('Fe', 0.2984)
barrel.add_element('Ni', 0.0613 )
barrel.add_element('Cr', 0.0608)
barrel.add_nuclide('Mn55', 0.0009)
barrel.add_element('Mo',0.0082)
barrel.temperature=shieldtemp
barrel.set_density('g/cm3', 3.92)


# In[20]:


# Geometry definitions

fuel_or = openmc.ZCylinder(r=0.3501) 
clad_ir = openmc.ZCylinder(r=0.3501) 
clad_or = openmc.ZCylinder(r=0.4057) 

top = openmc.ZPlane(z0=+0.5)
bottom = openmc.ZPlane(z0=-0.5) 

fuel_region = -fuel_or & -top & +bottom
clad_region = +clad_ir & -clad_or  & -top & +bottom
coolant_region = +clad_or & -top & +bottom
 
fuel_cell=openmc.Cell(fill=fuel, region=fuel_region)
clad_cell = openmc.Cell(fill=HT_9, region=clad_region)
sodium_cell = openmc.Cell(fill=coolant, region=coolant_region)

fuel_univ=openmc.Universe(cells=[fuel_cell, clad_cell, sodium_cell])

# Creating filling for emtpy space in the core

sodium_cool_cell = openmc.Cell(fill=coolant)
sodium_cool_u = openmc.Universe(cells=(sodium_cool_cell,))


# In[21]:


# Create a hexagonal fuel assembly
fuel_assembly = openmc.model.hexagonal_prism(edge_length=7.8976, orientation='x') 
fuel_assembly_cladding=openmc.model.hexagonal_prism(edge_length=8.2460, orientation='x')

# Creating universe for outer fuel pins
fuel_ass_inner_cells = openmc.Cell(fill=coolant, region=fuel_assembly & +bottom & -top)
fuel_ass_inner_univ = openmc.Universe(cells=[fuel_ass_inner_cells])

# Creating universe for assembly duct
fuel_ass_duct_cells = openmc.Cell(name='fuel_ass_duct_cells',fill=HT_9, region= ~fuel_assembly & fuel_assembly_cladding & +bottom & -top)
fuel_ass_duct_univ = openmc.Universe(cells=[fuel_ass_duct_cells])

# Define a lattice for fuel assemblies

fuel_ass_lattice = openmc.HexLattice(name='fuel_ass_lattice')
fuel_ass_lattice.orientation='x'
fuel_ass_lattice.center=(0,0)
fuel_ass_lattice.pitch=[0.9134]

# Create rings of fuel universes that will fill the lattice
inone = [fuel_univ]*48
intwo = [fuel_univ]*42
inthree = [fuel_univ]*36
infour = [fuel_univ]*30
infive = [fuel_univ]*24
insix = [fuel_univ]*18
inseven = [fuel_univ]*12
ineight = [fuel_univ]*6
innine = [fuel_univ]*1
fuel_ass_lattice.universes = \
[inone,
 intwo,
 inthree,
 infour,
 infive,
 insix,
 inseven,
 ineight,
 innine]

# Set the outer of the hexagonal lattice with sodium
fuel_ass_lattice.outer= fuel_ass_inner_univ

# Fill a cell with the lattice. This cell is filled with the lattice and contained within the prism.
fuel_ass_cells = openmc.Cell(fill=fuel_ass_lattice, region=fuel_assembly & -top & +bottom)
out_fuel_ass_cells  = openmc.Cell(fill=coolant, region=~fuel_assembly_cladding & -top & +bottom)


# Create a universe that contains both 
fuel_ass_univ = openmc.Universe(cells=[fuel_ass_cells, fuel_ass_duct_cells, out_fuel_ass_cells])


fuel_ass_univ.plot(origin = (0,0,0), pixels=(500, 500), width = (20.,20.), color_by = 'material')


# Create a hexagonal reflector assembly

reflector_assembly = openmc.model.hexagonal_prism(edge_length=7.8976, orientation='x')
reflector_assembly_cladding = openmc.model.hexagonal_prism(edge_length=8.2460, orientation='x')

#Create a reflector rod
ref=openmc.model.ZCylinder(r=0.7068)
ref_rod_cells=[openmc.Cell(name='ref_rod_cells', fill=reflector, region=-ref & -top & + bottom), openmc.Cell(fill=coolant, region=+ref & -top & +bottom)]
ref_rod_univ=openmc.Universe(name='ref_rod', cells=ref_rod_cells)

#Create the filling between the rods with coolant
reflector_ass_inner_cells=openmc.Cell(name='ref_rod_out_cells', fill=coolant, region=reflector_assembly & +bottom & -top)
reflector_ass_inner_univ= openmc.Universe(name='ref_rod_out', cells=[reflector_ass_inner_cells])

# Creating universe for assembly duct
ref_ass_duct_cells = openmc.Cell(name=' ref_duct_cells', fill=HT_9, region= ~reflector_assembly & reflector_assembly_cladding & +bottom & -top)
ref_ass_duct_univ = openmc.Universe(name='ref_ass_duct', cells=[ref_ass_duct_cells])

ref_ass_lattice=openmc.HexLattice(name='ref_ass_lattice')
ref_ass_lattice.orientation='x'
ref_ass_lattice.center=(0,0)
ref_ass_lattice.pitch=[1.4151]
ref_ass_lattice.universes=\
[
    [ref_rod_univ]*30,
    [ref_rod_univ]*24,
    [ref_rod_univ]*18,
    [ref_rod_univ]*12,
    [ref_rod_univ]*6,
    [ref_rod_univ]    
]

ref_ass_lattice.outer=reflector_ass_inner_univ
ref_ass_cells=openmc.Cell(name='ref_ass_cells', fill=ref_ass_lattice, region=reflector_assembly & +bottom & -top)
out_ref_ass_cells  = openmc.Cell(fill=coolant, region=~reflector_assembly_cladding & -top & +bottom)
ref_ass_univ=openmc.Universe(name='ref_ass_univ', cells=[ref_ass_cells, ref_ass_duct_cells, out_ref_ass_cells])


# Create a control assembly with central sodium bond

control_assembly=openmc.model.hexagonal_prism(edge_length=7.8976, orientation='x')
control_assembly_cladding=openmc.model.hexagonal_prism(edge_length=8.2460, orientation='x')
control_ass_cell=openmc.Cell(fill=coolant, region=control_assembly & -top & +bottom)
control_ass_clad_cell=openmc.Cell(fill=HT_9, region=~control_assembly & control_assembly_cladding & -top & +bottom)
out_control_ass_cells  = openmc.Cell(fill=coolant, region=~control_assembly & -top & +bottom)
control_univ = openmc.Universe(cells=[control_ass_cell,control_ass_clad_cell, out_control_ass_cells])


# Create a hexagonal shield assembly, in the homogeneous model

shield_assembly = openmc.model.hexagonal_prism(edge_length=8.2460, orientation='x')
shield_cell = openmc.Cell(fill=shield, region=shield_assembly & -top & +bottom)
out_shield_ass_cells  = openmc.Cell(fill=coolant, region=~shield_assembly & -top & +bottom)
shield_univ = openmc.Universe(cells=[shield_cell,out_shield_ass_cells])

# Create a hexagonal barrel assembly, in the homogeneous model

barrel_assembly = openmc.model.hexagonal_prism(edge_length=8.2460, orientation='x')
barrel_cell = openmc.Cell(fill=barrel, region=barrel_assembly & -top & +bottom)
barrel_univ = openmc.Universe(cells=[barrel_cell])


# In[26]:


# Creating the core

rad1=6*14.2826 #34.0 #cm
rad2= rad1+ 14.2826 # core barrel tickness is 14cm

core_c1=openmc.ZCylinder(r=rad1)
core_c2=openmc.ZCylinder(r=rad2, boundary_type='vacuum')

core_z0 = openmc.ZPlane(z0=-0.5,boundary_type='reflective')
core_z1 = openmc.ZPlane(z0=+0.5,boundary_type='reflective')

core1=-core_c1 & +core_z0 & -core_z1
core2=+core_c1 & -core_c2 & +core_z0 & -core_z1


core_inner_cells=openmc.Cell(fill=barrel, region=core1)
core_inner_univ=openmc.Universe(cells=[core_inner_cells])

absc=control_univ
fc=fuel_ass_univ
rc=ref_ass_univ
s=shield_univ
b=barrel_univ


# Define the core lattice

core_lattice = openmc.HexLattice(name='core_lattice')
core_lattice.center = (0., 0.)
core_lattice.orientation='y'
core_lattice.pitch = [14.6850]
core_lattice.outer = core_inner_univ



# In[15]:


# Create rings of fuel universes that will fill the lattice
shield_one = [s]*30
ref_one = [rc] * 24
ref_two = [rc]*18
fuel_one=[fc]*12
fuel_two=[fc]*6
abs_one=[absc]

core_lattice.universes = [shield_one, ref_one,ref_two,fuel_one,fuel_two, abs_one]

core_cells=openmc.Cell(fill=core_lattice, region=core1)
barrel_ring_cells=openmc.Cell(fill=barrel, region=core2)

core_univ=openmc.Universe(cells=[core_cells, barrel_ring_cells])


# Export geometry

geom = openmc.Geometry(core_univ)
geom.export_to_xml()
mats = list(geom.get_all_materials().values())
openmc.Materials(mats).export_to_xml()


# Settings
lower_left = np.array([-rad1, -rad1, -0.5])
upper_right = np.array([rad1, rad1, 0.5])


box = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
src = openmc.Source(space=box)
settings = openmc.Settings()
settings.source = src
settings.inactive = 20
settings.batches = 200
settings.particles = 5000
settings.temperature={'method': 'nearest','tolerance':1000.0}

settings.export_to_xml()

# Plot
ref_ass_univ.plot(origin = (0,0,0), pixels=(500, 500), width = (20.,20.), color_by = 'material')

control_univ.plot(origin = (0,0,0), pixels=(500, 500), width = (20.,20.), color_by = 'material')

core_univ.plot(origin = (0,0,0), pixels=(1000, 1000), width = (200,200), color_by = 'material')


########################################################################################################################
                                        # DEPLETION
########################################################################################################################

import openmc.deplete

# In 2D the height has to be set to 1.0

# New calculation for the kernel volume

N_holes= 217 # Number of fuel holes per fuel blocks
volfuel= 18*N_holes*((np.pi*0.3501**2*1.0))

print('Volume fuel :', volfuel)

fuel.volume = volfuel

Modelicius = openmc.model.Model(geometry=geom, materials=mats, settings=settings)
Modelicius.export_to_xml()

operator = openmc.deplete.Operator(Modelicius, chain_file, diff_burnable_mats=False)

fuel_block_edge=8.2460
power_density=68.8 #W/cm3
linear_power=250000 #W/cm
print('Linear power:' +str(linear_power))

# Calculating time steps for burn-up of 120 MWd/kgU
BU_final=120 #MWd/kgU
Max_step= 2*operator.heavy_metal/linear_power*1000
print('Max_step :', Max_step)
N_days=(BU_final*operator.heavy_metal*1000)/linear_power
end_step=N_days-9800

time_operation=[0.25]+[6.00]+[6.25]+[12.50]+[25.00]+[25.00]+[25.00]+([25.00]*16) + ([200]*10)+([365]*20)+([end_step])

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
power=[linear_power]*len(time_operation)
cooling=[0]*len(time_cooling)
power=power+cooling

print("Running time: "+str(N_days) + ' nombre de jours simulations '+str (sum(time_operation)))

integrator = openmc.deplete.CECMIntegrator(operator, time_steps, power, timestep_units='d')
openmc.run()
print("operator.heavy_metal : " + str(operator.heavy_metal))

start_time = time.time()
#integrator.integrate()
print("--- %s seconds ---" % (time.time() - start_time))



