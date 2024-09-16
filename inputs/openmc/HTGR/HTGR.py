import numpy as np
import openmc
import time


# Core parameters


fueltemp = 900 #K
lbptemp = 600 #K
modtemp = 600 #K
refltemp = 600 #K
bartemp = 600 #K
cooltemp = 600 #K

R = 8.314 #m3Pa /K/mol
T = cooltemp
M = 4 #g/mol
P = 6.39e6 #Pa
rho = M * P / (R*T) / 1e6 #g/cm3

fuelpf = 0.35
lbppf = 0.109





# Material

# Path to the cross section file
#openmc.Materials.cross_sections='/home/shared/endfb80_hdf5/cross_sections.xml'

kernel = openmc.Material(name='kernel')
kernel.set_density('g/cm3', 10.5)
kernel.add_nuclide('U235', 0.156672)
kernel.add_nuclide('U238', 0.843328)
kernel.add_nuclide('O16', 1.5)
kernel.add_element('C', 0.5)
kernel.temperature = fueltemp

buff= openmc.Material(name='buffer')
buff.set_density('g/cm3', 1.0)
buff.add_element('C', 1.0)
buff.add_s_alpha_beta('c_Graphite')
buff.temperature = fueltemp

iPyC = openmc.Material(name='iPyC')
iPyC.set_density('g/cm3', 1.9)
iPyC.add_element('C', 1.0)
iPyC.add_s_alpha_beta('c_Graphite')
iPyC.temperature = fueltemp

SiC = openmc.Material(name='SiC')
SiC.set_density('g/cm3', 3.2)
SiC.add_element('C', 0.5)
SiC.add_element('Si', 0.5)
SiC.temperature = fueltemp

oPyC = openmc.Material(name='oPyC')
oPyC.set_density('g/cm3', 1.9)
oPyC.add_element('C', 1.0)
oPyC.add_s_alpha_beta('c_Graphite')
oPyC.temperature = fueltemp

compact_matrix = openmc.Material(name='compact_matrix')
compact_matrix.set_density('g/cm3', 1.6479575)
compact_matrix.add_element('C',    0.999999540, percent_type='wo')
compact_matrix.add_nuclide('B10',  0.000000150, percent_type='wo')
compact_matrix.add_nuclide('Si28', 0.000000100, percent_type='wo')
compact_matrix.add_nuclide('Ca40', 0.000000080, percent_type='wo')
compact_matrix.add_nuclide('Fe56', 0.000000060, percent_type='wo')
compact_matrix.add_nuclide('Al27', 0.000000012, percent_type='wo')
compact_matrix.add_nuclide('K39',  0.000000040, percent_type='wo')
compact_matrix.add_nuclide('V51',  0.000000018, percent_type='wo')
compact_matrix.add_s_alpha_beta('c_Graphite')
compact_matrix.temperature = fueltemp

compact_graphite = openmc.Material(name='compact_graphite')
compact_graphite.set_density('g/cm3', 1.6479575)
compact_graphite.add_element('C',    0.999999540, percent_type='wo')
compact_graphite.add_nuclide('B10',  0.000000150, percent_type='wo')
compact_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
compact_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
compact_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
compact_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
compact_graphite.add_nuclide('K39',  0.000000040, percent_type='wo')
compact_graphite.add_nuclide('V51',  0.000000018, percent_type='wo')
compact_graphite.add_s_alpha_beta('c_Graphite')
compact_graphite.temperature = fueltemp

block_graphite = openmc.Material(name='block_graphite')
block_graphite.set_density('g/cm3', 1.85)
block_graphite.add_element('C',    0.999999540, percent_type='wo')
block_graphite.add_nuclide('B10',  0.000000150, percent_type='wo')
block_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
block_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
block_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
block_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
block_graphite.add_nuclide('K39',  0.000000040, percent_type='wo')
block_graphite.add_nuclide('V51',  0.000000018, percent_type='wo')
block_graphite.add_s_alpha_beta('c_Graphite')
block_graphite.temperature = modtemp

"""
B4C = openmc.Material(name='B4C')
B4C.set_density('g/cm3', 2.47)
B4C.add_element('B', 4.0)
B4C.add_element('C', 1.0)
B4C.add_s_alpha_beta('c_Graphite')
B4C.temperature = lbptemp
B4C.depletable = True

PyC = openmc.Material(name='PyC')
PyC.set_density('g/cm3', 1.87)
PyC.add_element('C', 1.0)
PyC.add_s_alpha_beta('c_Graphite')
PyC.temperature = lbptemp

buf = openmc.Material(name='buf')
buf.set_density('g/cm3', 1.0)
buf.add_element('C', 1.0)
buf.add_s_alpha_beta('c_Graphite')
buf.temperature = lbptemp

lbp_matrix = openmc.Material(name='lbp_matrix')
lbp_matrix.set_density('g/cm3', 1.37)
lbp_matrix.add_element('C',    0.999999540, percent_type='wo')
lbp_matrix.add_nuclide('B10',  0.000000150, percent_type='wo')
lbp_matrix.add_nuclide('Si28', 0.000000100, percent_type='wo')
lbp_matrix.add_nuclide('Ca40', 0.000000080, percent_type='wo')
lbp_matrix.add_nuclide('Fe56', 0.000000060, percent_type='wo')
lbp_matrix.add_nuclide('Al27', 0.000000012, percent_type='wo')
lbp_matrix.add_nuclide('K39',  0.000000040, percent_type='wo')
lbp_matrix.add_nuclide('V51',  0.000000018, percent_type='wo')
lbp_matrix.add_s_alpha_beta('c_Graphite')
lbp_matrix.temperature = lbptemp

lbp_graphite = openmc.Material(name='lbp_graphite')
lbp_graphite.set_density('g/cm3', 1.37)
lbp_graphite.add_element('C',    0.999999540, percent_type='wo')
lbp_graphite.add_nuclide('B10',  0.000000150, percent_type='wo')
lbp_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
lbp_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
lbp_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
lbp_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
lbp_graphite.add_nuclide('K39',  0.000000040, percent_type='wo')
lbp_graphite.add_nuclide('V51',  0.000000018, percent_type='wo')
lbp_graphite.add_s_alpha_beta('c_Graphite')
lbp_graphite.temperature = lbptemp
"""
coolant = openmc.Material(name='coolant')
coolant.set_density('g/cm3', rho) #0.01046)
coolant.add_nuclide('He4', 1.0)
coolant.temperature = cooltemp

reflector_graphite = openmc.Material(name='reflector_graphite')
reflector_graphite.set_density('g/cm3', 1.85)
reflector_graphite.add_element('C',    0.999999540, percent_type='wo')
reflector_graphite.add_nuclide('B10',  0.000000150, percent_type='wo')
reflector_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
reflector_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
reflector_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
reflector_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
reflector_graphite.add_nuclide('K39',  0.000000040, percent_type='wo')
reflector_graphite.add_nuclide('V51',  0.000000018, percent_type='wo')
reflector_graphite.add_s_alpha_beta('c_Graphite')
reflector_graphite.temperature = refltemp

permanent_graphite = openmc.Material(name='permanent_graphite')
permanent_graphite.set_density('g/cm3', 1.85)
permanent_graphite.add_element('C',    0.999999540, percent_type='wo')
permanent_graphite.add_nuclide('B10',  0.000000150, percent_type='wo')
permanent_graphite.add_nuclide('Si28', 0.000000100, percent_type='wo')
permanent_graphite.add_nuclide('Ca40', 0.000000080, percent_type='wo')
permanent_graphite.add_nuclide('Fe56', 0.000000060, percent_type='wo')
permanent_graphite.add_nuclide('Al27', 0.000000012, percent_type='wo')
permanent_graphite.add_nuclide('K39',  0.000000040, percent_type='wo')
permanent_graphite.add_nuclide('V51',  0.000000018, percent_type='wo')
permanent_graphite.add_s_alpha_beta('c_Graphite')
permanent_graphite.temperature = refltemp

barrel = openmc.Material(name='barrel')
barrel.set_density('g/cm3', 7.94)
barrel.add_element('Ni', 32.5)
barrel.add_element('Cr', 21)
barrel.add_element('Fe', 39.5)
barrel.add_element('Al', 0.3)
barrel.add_element('Ti', 0.3)
barrel.add_element('C', 0.08)
barrel.temperature = bartemp


##########################################################################################
# Creating homogenenous TRISO particles
# DEtermining the corresponding radius

hete_triso_radius = [2125e-5, 3125e-5, 3475e-5, 3825e-5, 4225e-5]
hete_lbp_radius = [1000e-5, 1180e-5, 1410e-5]

r0 = 0.5 #fix ratio of each layer equally except kernel
r1 = 0.583 #ratio of compact_matrix
r2 = 0.583 #ratio of lbp_matrix

num_triso = 7091
num_lbp = 52742

hete_kernel = np.pi*(4*hete_triso_radius[0]**3/3)*num_triso*15*10
hete_buffer = np.pi*(4*(hete_triso_radius[1]**3-hete_triso_radius[0]**3)/3)*num_triso*15*10
hete_iPyC = np.pi*(4*(hete_triso_radius[2]**3-hete_triso_radius[1]**3)/3)*num_triso*15*10
hete_SiC = np.pi*(4*(hete_triso_radius[3]**3-hete_triso_radius[2]**3)/3)*num_triso*15*10
hete_oPyC = np.pi*(4*(hete_triso_radius[4]**3-hete_triso_radius[3]**3)/3)*num_triso*15*10
hete_PyC = hete_iPyC + hete_oPyC 

hete_triso = hete_kernel + hete_buffer + hete_PyC + hete_SiC
hete_matrix = (np.pi*(0.635**2)*793) - hete_triso

triso_radius = []

R0  = ((hete_matrix*r1)/793/np.pi)**(0.5)
R1  = ((hete_matrix*r1 + hete_oPyC*r0)/793/np.pi)**0.5
R2  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0)/793/np.pi)**0.5
R3  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0)/793/np.pi)**0.5
R4  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0 + hete_buffer*r0)/793/np.pi)**0.5
R5  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0 + hete_buffer*r0 + hete_kernel)/793/np.pi)**0.5
R6  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC*r0 + hete_buffer    + hete_kernel)/793/np.pi)**0.5
R7  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC*r0 + hete_iPyC    + hete_buffer    + hete_kernel)/793/np.pi)**0.5
R8  = ((hete_matrix*r1 + hete_oPyC*r0 + hete_SiC    + hete_iPyC    + hete_buffer    + hete_kernel)/793/np.pi)**0.5
R9  = ((hete_matrix*r1 + hete_oPyC    + hete_SiC    + hete_iPyC    + hete_buffer    + hete_kernel)/793/np.pi)**0.5
R10 = ((hete_matrix    + hete_oPyC    + hete_SiC    + hete_iPyC    + hete_buffer    + hete_kernel)/793/np.pi)**0.5

triso_radius.append(R0)
triso_radius.append(R1)
triso_radius.append(R2)
triso_radius.append(R3)
triso_radius.append(R4)
triso_radius.append(R5)
triso_radius.append(R6)
triso_radius.append(R7)
triso_radius.append(R8)
triso_radius.append(R9)
triso_radius.append(R10)
print("\ntriso_radius :")
print(triso_radius)


########### Creating TRISO ##########

triso_spheres = [openmc.model.ZCylinder(r=r) for r in triso_radius]
triso_z0 = openmc.ZPlane(z0=-1/2) # 2D simulation to simplify the design
triso_z1 = openmc.ZPlane(z0=1/2)


old_triso_cells = []
old_triso_univ = []
for i in range(1):
    old_triso_cells.append(\
                           [openmc.Cell(fill=compact_matrix,  region=                    -triso_spheres[0] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=oPyC,            region=+triso_spheres[0] & -triso_spheres[1] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=SiC,             region=+triso_spheres[1] & -triso_spheres[2] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=iPyC,            region=+triso_spheres[2] & -triso_spheres[3] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=buff,            region=+triso_spheres[3] & -triso_spheres[4] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=kernel,          region=+triso_spheres[4] & -triso_spheres[5] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=buff,            region=+triso_spheres[5] & -triso_spheres[6] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=iPyC,            region=+triso_spheres[6] & -triso_spheres[7] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=SiC,             region=+triso_spheres[7] & -triso_spheres[8] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=oPyC,            region=+triso_spheres[8] & -triso_spheres[9] & +triso_z0 & -triso_z1),
                            openmc.Cell(fill=compact_matrix,  region=+triso_spheres[9] & -triso_spheres[10]& +triso_z0 & -triso_z1)])
    old_triso_univ.append(openmc.Universe(cells=old_triso_cells[i]))


triso_cells = []
triso_outer_cells = []
triso_univ = []
for i in range(1):
    triso_cells.append(openmc.Cell(fill=old_triso_univ[i], region=-triso_spheres[-1] & +triso_z0 & -triso_z1))
    triso_outer_cells.append(openmc.Cell(fill=block_graphite, region=+triso_spheres[-1] | -triso_z0 | +triso_z1))
    triso_univ.append(openmc.Universe(cells=[triso_cells[i],triso_outer_cells[i]]))


# In[28]:


# Creating Fuel Block

fuel_block = openmc.model.hexagonal_prism(edge_length=20.6756,orientation='x')
fuel_block_z0 = openmc.ZPlane(z0=-1/2)
fuel_block_z1 = openmc.ZPlane(z0=1/2)

small_coolant = openmc.model.ZCylinder(r=0.635)
small_coolant_cells = [openmc.Cell(fill=coolant, region=-small_coolant),
                      openmc.Cell(fill=block_graphite, region=+small_coolant)]
small_coolant_univ = openmc.Universe(cells=small_coolant_cells)

large_coolant = openmc.model.ZCylinder(r=0.794)
large_coolant_cells = [openmc.Cell(fill=coolant, region=-large_coolant),
                      openmc.Cell(fill=block_graphite, region=+large_coolant)]
large_coolant_univ = openmc.Universe(cells=large_coolant_cells)

hole = openmc.model.ZCylinder(r=0.635)
hole_cells = [openmc.Cell(fill=block_graphite, region=-hole),
             openmc.Cell(fill=block_graphite, region=+hole)]
hole_univ = openmc.Universe(cells=hole_cells)

a = triso_univ[0]
b = large_coolant_univ
c = small_coolant_univ
d = hole_univ


fuel_block_inner_cells = openmc.Cell(fill=block_graphite, region=fuel_block & +fuel_block_z0 & -fuel_block_z1)
fuel_block_inner_univ = openmc.Universe(cells=[fuel_block_inner_cells])

# Modifying fuel block lattice by removing lbp poison
fuel_block_lattice = openmc.HexLattice(name='fuel_block_lattice')
fuel_block_lattice.orientation='x'
fuel_block_lattice.center = (0, 0)
fuel_block_lattice.pitch = [1.8796]
fuel_block_lattice.universes = \
[
    [a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a], #11
    [b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a], #10
    [a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b], #9
    [a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[b]+[a], #8
    [b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a], #7
    [a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b]+[a]+[b]+[a]+[a]+[b], #6
    [a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a]+[a]+[a]+[b]+[a], #5
    [b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a]+[b]+[a]+[a], #4
    [a]+[c]+[a]+[c]+[a]+[c]+[a]+[c]+[a]+[c]+[a]+[c], #3
    [d]*6, #2
    [d], #1
]

fuel_block_lattice.outer = fuel_block_inner_univ

fuel_block_cells = openmc.Cell(name='fuel_block_cells',fill=fuel_block_lattice, region=fuel_block & +fuel_block_z0 & -fuel_block_z1)
fuel_block_univ = openmc.Universe(cells=[fuel_block_cells])
fuel_columns_univ = openmc.Universe(cells=[fuel_block_cells]) 


# Creating reflector column 

reflector = openmc.model.hexagonal_prism(edge_length=20.6756,orientation='x')

reflector_z0 = openmc.ZPlane(z0=-0.5)
reflector_z1 = openmc.ZPlane(z0=0.5)

reflector_cells = openmc.Cell(fill=reflector_graphite, region=reflector & +reflector_z0 & -reflector_z1)
reflector_univ = openmc.Universe(cells=[reflector_cells])


# Creating None Column

none = openmc.model.hexagonal_prism(edge_length=20.6756,orientation='x')
none_z0 = openmc.ZPlane(z0=-0.5)
none_z1 = openmc.ZPlane(z0=0.5)

none_cells = openmc.Cell(fill=permanent_graphite,region=none & +none_z0 & -none_z1)
none_univ = openmc.Universe(cells=[none_cells])

# Creating Core -> modyfying Cristina's core -> no lbp, no RSC holes

rad1 = 70.0 # radius of the active core
rad2= 110.0 # outer radius for the reflector
rad3 = rad2+(304.9-297.3) # outer vessel radius

# Reflective boundary conditions for an infinitely tall core
core_z0 = openmc.ZPlane(z0=-0.5,boundary_type='reflective') 
core_z1 = openmc.ZPlane(z0=0.5,boundary_type='reflective')

core_c1 = openmc.ZCylinder(r=rad1)
core_c2 = openmc.ZCylinder(r=rad2)
core_c3=openmc.ZCylinder(r=rad3, boundary_type='vacuum' )

core1 = -core_c1 & +core_z0 & -core_z1
core2 = +core_c1 & -core_c2 & +core_z0 & -core_z1
core1_inner_cells = openmc.Cell(fill=permanent_graphite, region=core1)
core1_inner_univ = openmc.Universe(cells=[core1_inner_cells])

fc = fuel_columns_univ
r = reflector_univ
n = none_univ

core_lattice = openmc.HexLattice(name='core_lattice')
core_lattice.orientation= 'y'
core_lattice.center = (0, 0)
core_lattice.pitch = [17.905594838485538*2]
core_lattice.universes = \
[
    [fc]*12, #3
    [fc]*6, #2
    [fc], #1
]
core_lattice.outer = core1_inner_univ

core_cells = openmc.Cell(fill=core_lattice, region=core1)
reflector_cells=openmc.Cell(fill=permanent_graphite, region=core2)
barrel_cells = openmc.Cell(fill=barrel, region=+core_c2 & -core_c3 & +core_z0 & -core_z1)

core_univ = openmc.Universe(cells=[core_cells, reflector_cells, barrel_cells])


# Export geometry


geom = openmc.Geometry(core_univ)
geom.export_to_xml()
mats = list(geom.get_all_materials().values())
openmc.Materials(mats).export_to_xml()


# Settings

lower_left = np.array([-70.0, -70.0, -0.5])
upper_right = np.array([70.0, 70.0, 0.5])

# The number of particles could be increased for a better accuracy 
box = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
src = openmc.Source(space=box)
settings = openmc.Settings()
settings.source = src
settings.inactive = 20
settings.batches = 200
settings.particles = 5000
settings.temperature={'method': 'interpolation','range':(293.15,923.15)}

settings.export_to_xml()

#########################################################################################################################
                                             # DEPLETION
#########################################################################################################################
import openmc.deplete


# New calculation for the kernel volume

fuel_block_edge=20.6756
N_holes= 216 # Number of fuel holes per fuel blocks

volkernel= N_holes*((np.pi*(triso_spheres[5].r**2-triso_spheres[4].r**2)*1.0))*16 -6*69*((np.pi*(triso_spheres[5].r**2-triso_spheres[4].r**2)*1.0))
print(volkernel)

kernel.volume = volkernel

Modelicius = openmc.model.Model(geometry=geom, materials=mats, settings=settings)
Modelicius.export_to_xml()

chain_file = '/home/chloec/Sources/chain_endfb71_pwr.xml'

operator = openmc.deplete.Operator(Modelicius, chain_file, diff_burnable_mats=False)

# Input power
power_density=5 #W/cm3
linear_power=power_density*(np.pi*(rad1**2)) #W/cm
print('Linear power:' +str(linear_power))

# Calculating time steps for burn-up of 60, 80, 1000 MWd/kgU
BU_final=80 #MWd/kgU # can be modified for other burnup
N_days=(BU_final*operator.heavy_metal*1000)/linear_power
end_step=N_days-1523

time_operation = ([3/24]*4)+([4/24]*3)+([6/24]*4)+([12/24]*2)+([1])+([3]*1)+([7]*3)+([14]*6)+([28]*8)+([29])+\
             ([31])+([31])+ ([31])+([30])+([31])+([30])+([31])+([31])+([30])+([31]) +\
             ([59])+([61])+([61])+([62])+([61])+([61])+\
             ([60])+([61])+ ([61])+([62])+([61])+([61])+\
             ([120])+([end_step])

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

print("operator.heavy_metal : " + str(operator.heavy_metal))


#############################################################################################
                                    #PLOTS
#############################################################################################
allcolors = {buff:'blue', iPyC:'yellow', SiC:'green', oPyC:'yellow', \
             compact_matrix:'cyan',compact_graphite:'purple', block_graphite:'yellow',\
              coolant:'lightblue', permanent_graphite:'green', \
              barrel:'black', kernel:'red', reflector_graphite:'grey',}

core_univ.plot(origin=(0.0, 0.0, 0.0), width=(250, 250), pixels=(1000,1000), basis='xy',color_by = 'material', colors=allcolors)

o = openmc.Plot.from_geometry(geom)
o.filename = 'core000'
o.origin = [0,0,0]
o.width = [250,250]
o.pixels = (10500,10500)
o.color_by = 'material'
o.basis = 'xy'
o.colors = allcolors



w = openmc.Plot.from_geometry(geom)
w.filename = 'core008'
w.origin = [0, 17.905594838485538*2,0]
w.width = [17.905594838485538*2+1,17.905594838485538*2+1]
w.pixels = (8000,8000)
w.color_by = 'material'
w.basis = 'xy'
w.colors = allcolors

plot_file = openmc.Plots((o,w))

plot_file.export_to_xml()
openmc.plot_geometry()

##################################################################################################
#                            DEPLETION
##################################################################################################

start_time = time.time()
integrator.integrate()
print("--- %s seconds ---" % (time.time() - start_time))

