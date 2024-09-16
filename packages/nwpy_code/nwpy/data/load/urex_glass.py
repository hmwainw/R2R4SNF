################################################################################
"""
    Last modified on May 16, 2019
        
    @author: matz
        
    comment: data file containing waste form loading information for glass 
             waste from urex aqeuous reprocessing, along with stoichiometry
             for the oxides of waste nuclides.
    
"""
################################################################################
from __future__ import absolute_import
import numpy as np
################################################################################
# GLASS from aqueous reprocessing (UREX/+)
canister = {}
canister['Diameter'] = 0.6096 # m
canister['Length'] = 4.572 # m
canister['Mass limit'] = 2.9e6 # grams
canister['Thickness'] = 0.06 # 1cm for canister, 5cm overpack thickness
################################################################################
# Definitions:
#   Mw = Mass of waste
#   Mwo = Mass of waste-oxides per canister
#   Mg = Mass of glass frit (waste form matrix/diluent) per canister
#   xo_(species) = Mass fraction of species in waste form (basis waste-oxides)
#   Ho = Decay heat per unit mass waste-oxide (i.e. power/gram)
#-------------------------------------------------------------------------------
# Constraints (6-8 require input from stream data)
# 1. Nonzero waste mass         Mwo > 0
# 2. Nonzero glass frit mass    Mg > 0
# 3. Mass limit: 2900 kg        Mwo + Mg <= 2900 kg (2.9e6 g)
# 4. Maximum loading: 0.45      Mwo / (Mwo + Mg) <= 0.45
# 5. Heat: 14000 kW/canister    Ho * Mwo <= 14 kW/canister
# 6. MoO3 fraction: 0.025       xo_(MoO3) * Mwo <= 0.025 * (Mwo + Mg)
# 7. Noble metal (RE) fraction  xo_(NM) * Mwo <= 0.03 * (Mwo + Mg)
# X. Minimum loading: 0.15      Mwo / (Mwo + Mg) >= 0.15
################################################################################
# CONSTRAINTS
glass = {}
glass['ineq'] = np.zeros((7, 3))
glass['ineq'][0] = [-1.0,     0.0,   0.0]
glass['ineq'][1] = [0.0,      -1.0,  0.0]
glass['ineq'][2] = [1.0,      1.0,   canister['Mass limit']]
glass['ineq'][3] = [1.0-0.45, -0.45, 0.0]
glass['ineq'][4] = [0.0,      0.0,   14000.0]
glass['ineq'][5] = [0.0,      0.0,   0.025]
glass['ineq'][6] = [0.0,      0.0,   0.03]
# glass['ineq'][3] = [0.15-1.0, 0.15,  0.0]
################################################################################
# KEY
glass['key'] = []
glass['key'].append('Nonzero waste mass')
glass['key'].append('Nonzero frit mass')
glass['key'].append('Max glass mass')
#glass['key'].append('Min waste loading')
glass['key'].append('Max waste loading')
glass['key'].append('Max heat')
glass['key'].append('Max MoO3 fraction')
glass['key'].append('Max rare earth fraction')
################################################################################
# OXIDE STOICHIOMETRY
################################################################################
# NOTES:
# Source: 2006 Ahn linear programming paper
# Moble metals: Rh, Ru, Ag, Pd
# Not in glass (therefore not here):
#	- halogens: F(9), Cl(17), Br(35), I(53), At(85)
#	- noble gases: He(2), Ne(10), Ar(18), Kr(36), Xe(54), Rn(86)
#	- other gases or volatile elements: H(1), C(6), N(7)
# Adjacent to the dictionary...
# * indicates oxide stoichiometry inferred by periodic table column
# ** indicates "educated" guess
#-------------------------------------------------------------------------------
oxide = {}
oxide['li'] = [2, 1]	# 3 Li2O
oxide['be'] = [1, 1]	# 4 BeO
oxide['b'] = [2, 3]		# 5 B2O3
oxide['c'] = [1, 2]		# 6 CO2
oxide['na'] = [2, 1]	# 11 Na2O
oxide['mg'] = [1, 1]	# 12 MgO
oxide['al'] = [2, 3]	# 13 Al2O3
oxide['si'] = [1, 2]	# 14 SiO2
oxide['p'] = [2, 5]		# 15 P2O5
oxide['s'] = [1, 2]		# 16 SO2	*
oxide['k'] = [2, 1]		# 19 K2O
oxide['ca'] = [1, 1]	# 20 CaO
oxide['sc'] = [2, 3]	# 21 Sc2O3
oxide['ti'] = [1, 2]	# 22 TiO2	*
oxide['v'] = [2, 5]		# 23 V2O5	*
oxide['cr'] = [2, 3]	# 24 Cr2O3
oxide['mn'] = [1, 2]	# 25 MnO2	*
oxide['fe'] = [2, 3]	# 26 Fe2O3
oxide['co'] = [2, 3]	# 27 Co2O3	*
oxide['ni'] = [1, 1]	# 28 NiO	*
oxide['cu'] = [2, 1]	# 29 Cu2O	*
oxide['zn'] = [1, 1]	# 30 ZnO	*
oxide['ga'] = [2, 3]	# 31 Ga2O3	*
oxide['ge'] = [1, 1]	# 32 GeO
oxide['as'] = [2, 3]	# 33 As2O3
oxide['se'] = [1, 2]	# 34 SeO2
oxide['rb'] = [2, 1]	# 37 Rb2O
oxide['sr'] = [1, 1]	# 38 SrO
oxide['y'] = [2, 3]		# 39 Y2O3
oxide['zr'] = [1, 2]	# 40 ZrO2
oxide['nb'] = [2, 5]	# 41 Nb2O5
oxide['mo'] = [1, 3]	# 42 MoO3
oxide['tc'] = [1, 2]	# 43 TcO2
oxide['ru'] = [1, 2]	# 44 RuO2
oxide['rh'] = [2, 3]	# 45 Rh2O3
oxide['pd'] = [1, 1]	# 46 PdO
oxide['ag'] = [2, 1]	# 47 Ag2O
oxide['cd'] = [1, 1]	# 48 CdO
oxide['in'] = [2, 3]	# 49 In2O3
oxide['sn'] = [1, 2]	# 50 SnO2
oxide['sb'] = [2, 3]	# 51 Sb2O3
oxide['te'] = [1, 2]	# 52 TeO2
oxide['cs'] = [2, 1]	# 55 Cs2O
oxide['ba'] = [1, 1]	# 56 BaO
oxide['la'] = [2, 3]	# 57 La2O3
oxide['ce'] = [1, 2]	# 58 CeO2
oxide['pr'] = [6, 11]	# 59 Pr6O11
oxide['nd'] = [2, 3]	# 60 Nd2O3
oxide['pm'] = [2, 3]	# 61 Pm2O3
oxide['sm'] = [2, 3]	# 62 Sm2O3
oxide['eu'] = [2, 3]	# 63 Eu2O3
oxide['gd'] = [2, 3]	# 64 Gd2O3
oxide['tb'] = [2, 3]	# 65 Tb2O3
oxide['dy'] = [2, 3]	# 66 Dy2O3
oxide['ho'] = [2, 3]	# 67 Ho2O3
oxide['er'] = [2, 3]	# 68 Er2O3	**
oxide['tm'] = [2, 3]	# 69 Tm2O3	**
oxide['yb'] = [2, 3]	# 70 Yb2O3	**
oxide['lu'] = [2, 3]	# 71 Lu2O3	**
oxide['hf'] = [1, 2]	# 72 HfO2	*
oxide['ta'] = [2, 5]	# 73 Ta2O5	*
oxide['w'] = [1, 2]		# 74 WO2	**
oxide['re'] = [1, 4]	# 75 ReO4	**
oxide['os'] = [1, 2]	# 76 OsO2	**
oxide['ir'] = [1, 2]	# 77 IrO2	**
oxide['pt'] = [1, 1]	# 78 PtO	*
oxide['au'] = [2, 3]	# 79 Au2O3	**
oxide['hg'] = [1, 1]	# 80 HgO	*
oxide['tl'] = [2, 3]	# 81 Tl2O3	*
oxide['pb'] = [1, 1]	# 82 PbO	**
oxide['bi'] = [2, 3]	# 83 Bi2O3	*
oxide['po'] = [1, 2]	# 84 PoO2	*
oxide['fr'] = [2, 1]	# 87 Fr2O	*
oxide['ra'] = [1, 1]	# 88 RaO	*
oxide['ac'] = [2, 3]	# 89 Ac2O3	*
oxide['th'] = [1, 2]	# 90 ThO2	**
oxide['pa'] = [6, 11]	# 91 Pa6O11	*
oxide['u'] = [3, 8]		# 92 U3O8
oxide['np'] = [1, 2]	# 93 NpO2
oxide['pu'] = [1, 2]	# 94 PuO2
oxide['am'] = [2, 3]	# 95 Am2O3
oxide['cm'] = [2, 3]	# 96 Cm2O3
oxide['bk'] = [2, 3]	# 97 Bk2O3	*
oxide['cf'] = [2, 3]	# 98 Cf2O3	*
oxide['es'] = [2, 3]	# 99 Es2O3	*
################################################################################
