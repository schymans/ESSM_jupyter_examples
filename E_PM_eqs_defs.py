from essm.variables._core import BaseVariable, Variable
from essm.equations import Equation
from sympy import Abs, exp, Float, log, nsolve, Eq, sqrt, solve, Symbol
from sympy.physics.units import second, pascal, mole, joule, meter, watt, kelvin, kilogram
alpha_a = type('alpha_a', (Variable,), {'__doc__': """Thermal diffusivity of dry air.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'\alpha_a', 'default': None, 'expr': None})
c_pa = type('c_pa', (Variable,), {'__doc__': """Specific heat of dry air.""", 'unit': joule/(kelvin*kilogram), 'assumptions': {'real': True},         'latex_name': r'c_{pa}', 'default': 1010.0, 'expr': None})
c_pamol = type('c_pamol', (Variable,), {'__doc__': """Molar specific heat of dry air.

    https://en.wikipedia.org/wiki/Heat_capacity#Specific_heat_capacity
    """, 'unit': joule/(kelvin*mole), 'assumptions': {'real': True},         'latex_name': r'c_{pa,mol}', 'default': 29.19, 'expr': None})
c_pv = type('c_pv', (Variable,), {'__doc__': """Specific heat of water vapour at 300 K.

    http://www.engineeringtoolbox.com/water-vapor-d_979.html
    """, 'unit': joule/(kelvin*kilogram), 'assumptions': {'real': True},         'latex_name': r'c_{pv}', 'default': 1864, 'expr': None})
C_wa = type('C_wa', (Variable,), {'__doc__': """Concentration of water in air.""", 'unit': mole/meter**3, 'assumptions': {'real': True},         'latex_name': r'C_{wa}', 'default': None, 'expr': None})
D_va = type('D_va', (Variable,), {'__doc__': """Binary diffusion coefficient of water vapour in air.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'D_{va}', 'default': None, 'expr': None})
g = type('g', (Variable,), {'__doc__': """Gravitational acceleration.""", 'unit': meter/second**2, 'assumptions': {'real': True},         'latex_name': r'g', 'default': 9.81, 'expr': None})
Gr = type('Gr', (Variable,), {'__doc__': """Grashof number.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Gr_L}', 'default': None, 'expr': None})
h_c = type('h_c', (Variable,), {'__doc__': """Average 1-sided convective heat transfer coefficient.""", 'unit': joule/(kelvin*meter**2*second), 'assumptions': {'real': True},         'latex_name': r'h_c', 'default': None, 'expr': None})
k_a = type('k_a', (Variable,), {'__doc__': """Thermal conductivity of dry air.""", 'unit': joule/(kelvin*meter*second), 'assumptions': {'real': True},         'latex_name': r'k_a', 'default': None, 'expr': None})
lambda_E = type('lambda_E', (Variable,), {'__doc__': """Latent heat of evaporation.""", 'unit': joule/kilogram, 'assumptions': {'real': True},         'latex_name': r'\lambda_E', 'default': 2450000.0, 'expr': None})
Le = type('Le', (Variable,), {'__doc__': """Lewis number.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Le}', 'default': None, 'expr': None})
M_air = type('M_air', (Variable,), {'__doc__': """Molar mass of air.

    http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
    """, 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_{air}', 'default': 0.02897, 'expr': None})
M_N2 = type('M_N2', (Variable,), {'__doc__': """Molar mass of nitrogen.""", 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_{N_2}', 'default': 0.028, 'expr': None})
M_O2 = type('M_O2', (Variable,), {'__doc__': """Molar mass of oxygen.""", 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_{O_2}', 'default': 0.032, 'expr': None})
M_w = type('M_w', (Variable,), {'__doc__': """Molar mass of water.""", 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_w', 'default': 0.018, 'expr': None})
nu_a = type('nu_a', (Variable,), {'__doc__': """Kinematic viscosity of dry air.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'\nu_a', 'default': None, 'expr': None})
Nu = type('Nu', (Variable,), {'__doc__': """Average Nusselt number over given length.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Nu_L}', 'default': None, 'expr': None})
P_a = type('P_a', (Variable,), {'__doc__': """Air pressure.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_a', 'default': None, 'expr': None})
Pr = type('Pr', (Variable,), {'__doc__': """Prandtl number (0.71 for air).""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Pr}', 'default': None, 'expr': None})
P_N2 = type('P_N2', (Variable,), {'__doc__': """Partial pressure of nitrogen.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{N2}', 'default': None, 'expr': None})
P_O2 = type('P_O2', (Variable,), {'__doc__': """Partial pressure of oxygen.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{O2}', 'default': None, 'expr': None})
P_wa = type('P_wa', (Variable,), {'__doc__': """Water vapour pressure in the atmosphere.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{wa}', 'default': None, 'expr': None})
P_was = type('P_was', (Variable,), {'__doc__': """Saturation water vapour pressure at air temperature.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{was}', 'default': None, 'expr': None})
R_d = type('R_d', (Variable,), {'__doc__': """Downwelling global radiation.""", 'unit': watt/meter**2, 'assumptions': {'real': True},         'latex_name': r'R_d', 'default': None, 'expr': None})
Re_c = type('Re_c', (Variable,), {'__doc__': """Critical Reynolds number for the onset of turbulence.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Re_c}', 'default': None, 'expr': None})
Re = type('Re', (Variable,), {'__doc__': """Average Reynolds number over given length.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Re_L}', 'default': None, 'expr': None})
rho_a = type('rho_a', (Variable,), {'__doc__': """Density of dry air.""", 'unit': kilogram/meter**3, 'assumptions': {'real': True},         'latex_name': r'\rho_a', 'default': None, 'expr': None})
R_u = type('R_u', (Variable,), {'__doc__': """Upwelling global radiation.""", 'unit': watt/meter**2, 'assumptions': {'real': True},         'latex_name': r'R_u', 'default': None, 'expr': None})
R_mol = type('R_mol', (Variable,), {'__doc__': """Molar gas constant.""", 'unit': joule/(kelvin*mole), 'assumptions': {'real': True},         'latex_name': r'R_{mol}', 'default': 8.314472, 'expr': None})
R_s = type('R_s', (Variable,), {'__doc__': """Solar shortwave flux per area.""", 'unit': joule/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'R_s', 'default': None, 'expr': None})
sigm = type('sigm', (Variable,), {'__doc__': """Stefan-Boltzmann constant.""", 'unit': joule/(kelvin**4*meter**2*second), 'assumptions': {'real': True},         'latex_name': r'\sigma', 'default': 5.67e-08, 'expr': None})
T0 = type('T0', (Variable,), {'__doc__': """Freezing point in Kelvin.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_0', 'default': 273.15, 'expr': None})
T_a = type('T_a', (Variable,), {'__doc__': """Air temperature.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_a', 'default': None, 'expr': None})
v_w = type('v_w', (Variable,), {'__doc__': """Wind velocity.""", 'unit': meter/second, 'assumptions': {'real': True},         'latex_name': r'v_w', 'default': None, 'expr': None})
x_N2 = type('x_N2', (Variable,), {'__doc__': """Mole fraction of nitrogen in dry air.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'x_{N2}', 'default': 0.79, 'expr': None})
x_O2 = type('x_O2', (Variable,), {'__doc__': """Mole fraction of oxygen in dry air.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'x_{O2}', 'default': 0.21, 'expr': None})
p_Dva1 = type('p_Dva1', (Variable,), {'__doc__': """Internal parameter of eq_Dva.""", 'unit': meter**2/(kelvin*second), 'assumptions': {'real': True},         'latex_name': r'p_1', 'default': 1.49e-07, 'expr': None})
p_Dva2 = type('p_Dva2', (Variable,), {'__doc__': """Internal parameter of eq_Dva.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'p_2', 'default': 1.96e-05, 'expr': None})
p_alpha1 = type('p_alpha1', (Variable,), {'__doc__': """Internal parameter of eq_alphaa.""", 'unit': meter**2/(kelvin*second), 'assumptions': {'real': True},         'latex_name': r'p_1', 'default': 1.32e-07, 'expr': None})
p_alpha2 = type('p_alpha2', (Variable,), {'__doc__': """Internal parameter of eq_alphaa.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'p_2', 'default': 1.73e-05, 'expr': None})
p_ka1 = type('p_ka1', (Variable,), {'__doc__': """Internal parameter of eq_ka.""", 'unit': joule/(kelvin**2*meter*second), 'assumptions': {'real': True},         'latex_name': r'p_1', 'default': 6.84e-05, 'expr': None})
p_ka2 = type('p_ka2', (Variable,), {'__doc__': """Internal parameter of eq_ka.""", 'unit': joule/(kelvin*meter*second), 'assumptions': {'real': True},         'latex_name': r'p_2', 'default': 0.00563, 'expr': None})
p_nua1 = type('p_nua1', (Variable,), {'__doc__': """Internal parameter of eq_nua.""", 'unit': meter**2/(kelvin*second), 'assumptions': {'real': True},         'latex_name': r'p_1', 'default': 9e-08, 'expr': None})
p_nua2 = type('p_nua2', (Variable,), {'__doc__': """Internal parameter of eq_nua.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'p_2', 'default': 1.13e-05, 'expr': None})
alpha_l = type('alpha_l', (Variable,), {'__doc__': """Leaf albedo, fraction of shortwave radiation reflected by the leaf.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'\alpha_l', 'default': None, 'expr': None})
a_s = type('a_s', (Variable,), {'__doc__': """Fraction of one-sided leaf area covered by stomata.

    (1 if stomata are on one side only, 2 if they are on both sides).
    """, 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'a_s', 'default': None, 'expr': None})
a_sh = type('a_sh', (Variable,), {'__doc__': """Fraction of projected area exchanging sensible heat with the air.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'a_{sh}', 'default': 2.0, 'expr': None})
C_wl = type('C_wl', (Variable,), {'__doc__': """Concentration of water in the leaf air space.""", 'unit': mole/meter**3, 'assumptions': {'real': True},         'latex_name': r'C_{wl}', 'default': None, 'expr': None})
E_lmol = type('E_lmol', (Variable,), {'__doc__': """Transpiration rate in molar units.""", 'unit': mole/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'E_{l,mol}', 'default': None, 'expr': None})
E_l = type('E_l', (Variable,), {'__doc__': """Latent heat flux from leaf.""", 'unit': joule/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'E_l', 'default': None, 'expr': None})
epsilon_l = type('epsilon_l', (Variable,), {'__doc__': """Longwave emmissivity of the leaf surface.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'\epsilon_l', 'default': 1.0, 'expr': None})
g_bw = type('g_bw', (Variable,), {'__doc__': """Boundary layer conductance to water vapour.""", 'unit': meter/second, 'assumptions': {'real': True},         'latex_name': r'g_{bw}', 'default': None, 'expr': None})
g_bwmol = type('g_bwmol', (Variable,), {'__doc__': """Boundary layer conductance to water vapour.""", 'unit': mole/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'g_{bw,mol}', 'default': None, 'expr': None})
g_sw = type('g_sw', (Variable,), {'__doc__': """Stomatal conductance to water vapour.""", 'unit': meter/second, 'assumptions': {'real': True},         'latex_name': r'g_{sw}', 'default': None, 'expr': None})
g_swmol = type('g_swmol', (Variable,), {'__doc__': """Stomatal conductance to water vapour.""", 'unit': mole/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'g_{sw,mol}', 'default': None, 'expr': None})
g_tw = type('g_tw', (Variable,), {'__doc__': """Total leaf conductance to water vapour.""", 'unit': meter/second, 'assumptions': {'real': True},         'latex_name': r'g_{tw}', 'default': None, 'expr': None})
g_twmol = type('g_twmol', (Variable,), {'__doc__': """Total leaf layer conductance to water vapour.""", 'unit': mole/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'g_{tw,mol}', 'default': None, 'expr': None})
H_l = type('H_l', (Variable,), {'__doc__': """Sensible heat flux from leaf.""", 'unit': joule/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'H_l', 'default': None, 'expr': None})
L_A = type('L_A', (Variable,), {'__doc__': """Leaf area.""", 'unit': meter**2, 'assumptions': {'real': True},         'latex_name': r'L_A', 'default': None, 'expr': None})
L_l = type('L_l', (Variable,), {'__doc__': """Leaf width as characteristic length scale for convection.""", 'unit': meter, 'assumptions': {'real': True},         'latex_name': r'L_l', 'default': None, 'expr': None})
P_wl = type('P_wl', (Variable,), {'__doc__': """Water vapour pressure inside the leaf.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{wl}', 'default': None, 'expr': None})
r_bw = type('r_bw', (Variable,), {'__doc__': """Boundary layer resistance to water vapour, inverse of $g_{bw}$.""", 'unit': second/meter, 'assumptions': {'real': True},         'latex_name': r'r_{bw}', 'default': None, 'expr': None})
r_sw = type('r_sw', (Variable,), {'__doc__': """Stomatal resistance to water vapour, inverse of $g_{sw}$.""", 'unit': second/meter, 'assumptions': {'real': True},         'latex_name': r'r_{sw}', 'default': None, 'expr': None})
r_tw = type('r_tw', (Variable,), {'__doc__': """Total leaf resistance to water vapour, $r_{bv} + r_{sv}$.""", 'unit': second/meter, 'assumptions': {'real': True},         'latex_name': r'r_{tw}', 'default': None, 'expr': None})
rho_al = type('rho_al', (Variable,), {'__doc__': """Density of air at the leaf surface.""", 'unit': kilogram/meter**3, 'assumptions': {'real': True},         'latex_name': r'\rho_{al}', 'default': None, 'expr': None})
R_la = type('R_la', (Variable,), {'__doc__': """Longwave radiation absorbed by leaf.""", 'unit': watt/meter**2, 'assumptions': {'real': True},         'latex_name': r'R_{la}', 'default': None, 'expr': None})
R_ll = type('R_ll', (Variable,), {'__doc__': """Longwave radiation away from leaf.""", 'unit': watt/meter**2, 'assumptions': {'real': True},         'latex_name': r'R_{ll}', 'default': None, 'expr': None})
R_ld = type('R_ld', (Variable,), {'__doc__': """Downwards emitted/reflected global radiation from leaf.""", 'unit': watt/meter**2, 'assumptions': {'real': True},         'latex_name': r'R_{ld}', 'default': None, 'expr': None})
R_lu = type('R_lu', (Variable,), {'__doc__': """Upwards emitted/reflected global radiation from leaf.""", 'unit': watt/meter**2, 'assumptions': {'real': True},         'latex_name': r'R_{lu}', 'default': None, 'expr': None})
T_l = type('T_l', (Variable,), {'__doc__': """Leaf temperature.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_l', 'default': None, 'expr': None})
T_w = type('T_w', (Variable,), {'__doc__': """Radiative temperature of objects surrounding the leaf.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_w', 'default': None, 'expr': None})
p_CC1 = type('p_CC1', (Variable,), {'__doc__': """Internal parameter of eq_Pwl.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'p_1', 'default': 611.0, 'expr': None})
p_CC2 = type('p_CC2', (Variable,), {'__doc__': """Internal parameter of eq_Pwl.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'p_2', 'default': 273.0, 'expr': None})
beta_B = type('beta_B', (Variable,), {'__doc__': """Bowen ratio (sensible/latent heat flux)""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'\beta_B', 'default': None, 'expr': None})
Delta_eTa = type('Delta_eTa', (Variable,), {'__doc__': """Slope of saturation vapour pressure at air temperature""", 'unit': pascal/kelvin, 'assumptions': {'real': True},         'latex_name': r'\Delta_{e_{Ta}}', 'default': None, 'expr': None})
epsilon = type('epsilon', (Variable,), {'__doc__': """Water to air molecular weight ratio (0.622)""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'\epsilon', 'default': 0.622, 'expr': None})
E_w = type('E_w', (Variable,), {'__doc__': """Latent heat flux from a wet surface""", 'unit': joule/(meter**2*second), 'assumptions': {'real': True},         'latex_name': r'E_w', 'default': None, 'expr': None})
f_u = type('f_u', (Variable,), {'__doc__': """Wind function in Penman approach, f(u) adapted to energetic units""", 'unit': joule/(meter**2*pascal*second), 'assumptions': {'real': True},         'latex_name': r'f_u', 'default': None, 'expr': None})
gamma_v = type('gamma_v', (Variable,), {'__doc__': """Psychrometric constant""", 'unit': pascal/kelvin, 'assumptions': {'real': True},         'latex_name': r'\gamma_v', 'default': None, 'expr': None})
n_MU = type('n_MU', (Variable,), {'__doc__': """n=2 for hypostomatous, n=1 for amphistomatous leaves""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'n_{MU}', 'default': None, 'expr': a_sh/a_s})
r_a = type('r_a', (Variable,), {'__doc__': """One-sided boundary layer resistance to heat transfer ($r_H$ in \citet[][P. 231]{monteith_principles_2013})""", 'unit': second/meter, 'assumptions': {'real': True},         'latex_name': r'r_a', 'default': None, 'expr': None})
r_v = type('r_v', (Variable,), {'__doc__': """One-sided leaf BL resistance to water vapour, \citep[][Eqs. 11.3 and 13.16]{monteith_principles_2013}""", 'unit': second/meter, 'assumptions': {'real': True},         'latex_name': r'r_{v}', 'default': None, 'expr': None})
r_s = type('r_s', (Variable,), {'__doc__': """Stomatal resistance to water vapour \citep[][P. 231]{monteith_principles_2013}""", 'unit': second/meter, 'assumptions': {'real': True},         'latex_name': r'r_s', 'default': None, 'expr': None})
S = type('S', (Variable,), {'__doc__': """Factor representing stomatal resistance in \citet{penman_physical_1952}""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'S', 'default': None, 'expr': None})
c_E = type('c_E', (Variable,), {'__doc__': """Latent heat transfer coefficient""", 'unit': joule/(meter**2*pascal*second), 'assumptions': {'real': True},         'latex_name': r'c_E', 'default': None, 'expr': None})
c_H = type('c_H', (Variable,), {'__doc__': """Sensible heat transfer coefficient""", 'unit': joule/(kelvin*meter**2*second), 'assumptions': {'real': True},         'latex_name': r'c_H', 'default': None, 'expr': None})
eq_Le = type('eq_Le', (Equation,), {'__doc__': """Le as function of alpha_a and D_va.

    (Eq. B3 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(Le, alpha_a/D_va)})
eq_Cwa = type('eq_Cwa', (Equation,), {'__doc__': """C_wa as a function of P_wa and T_a.

    (Eq. B9 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(C_wa, P_wa/(R_mol*T_a))})
eq_Nu_forced_all = type('eq_Nu_forced_all', (Equation,), {'__doc__': """Nu as function of Re and Re_c under forced conditions.

    (Eqs. B13--B15 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(Nu, -Pr**(1/3)*(-37*Re**(4/5) + 37*(Re + Re_c - Abs(Re - Re_c)/2)**(4/5) - 664*sqrt(Re + Re_c - Abs(Re - Re_c)/2))/1000)})
eq_Dva = type('eq_Dva', (Equation,), {'__doc__': """D_va as a function of air temperature.

    (Table A.3 in :cite:`monteith_principles_2007`)
    """, 'expr': Eq(D_va, T_a*p_Dva1 - p_Dva2)})
eq_alphaa = type('eq_alphaa', (Equation,), {'__doc__': """alpha_a as a function of air temperature.

    (Table A.3 in :cite:`monteith_principles_2007`)
    """, 'expr': Eq(alpha_a, T_a*p_alpha1 - p_alpha2)})
eq_ka = type('eq_ka', (Equation,), {'__doc__': """k_a as a function of air temperature.

    (Table A.3 in :cite:`monteith_principles_2007`)
    """, 'expr': Eq(k_a, T_a*p_ka1 + p_ka2)})
eq_nua = type('eq_nua', (Equation,), {'__doc__': """nu_a as a function of air temperature.

    (Table A.3 in :cite:`monteith_principles_2007`)
    """, 'expr': Eq(nu_a, T_a*p_nua1 - p_nua2)})
eq_rhoa_Pwa_Ta = type('eq_rhoa_Pwa_Ta', (Equation,), {'__doc__': """rho_a as a function of P_wa and T_a.

    (Eq. B20 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(rho_a, (M_N2*P_N2 + M_O2*P_O2 + M_w*P_wa)/(R_mol*T_a))})
eq_Pa = type('eq_Pa', (Equation,), {'__doc__': """Calculate air pressure.

    From partial pressures of N2, O2 and H2O, following Dalton's law of
    partial pressures.
    """, 'expr': Eq(P_a, P_N2 + P_O2 + P_wa)})
eq_PN2_PO2 = type('eq_PN2_PO2', (Equation,), {'__doc__': """Calculate P_N2 as a function of P_O2.

    It follows Dalton's law of partial pressures.
    """, 'expr': Eq(P_N2, P_O2*x_N2/x_O2)})
eq_PO2 = type('eq_PO2', (Equation,), {'__doc__': """Calculate P_O2 as a function of P_a, P_N2 and P_wa.""", 'expr': Eq(P_O2, (P_a*x_O2 - P_wa*x_O2)/(x_N2 + x_O2))})
eq_PN2 = type('eq_PN2', (Equation,), {'__doc__': """Calculate P_N2 as a function of P_a, P_O2 and P_wa.""", 'expr': Eq(P_N2, (P_a*x_N2 - P_wa*x_N2)/(x_N2 + x_O2))})
eq_rhoa = type('eq_rhoa', (Equation,), {'__doc__': """Calculate rho_a from T_a, P_a and P_wa.""", 'expr': Eq(rho_a, (x_N2*(M_N2*P_a - P_wa*(M_N2 - M_w)) + x_O2*(M_O2*P_a - P_wa*(M_O2 - M_w)))/(R_mol*T_a*x_N2 + R_mol*T_a*x_O2))})
eq_Rs_enbal = type('eq_Rs_enbal', (Equation,), {'__doc__': """Calculate R_s from energy balance.

    (Eq. 1 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(R_s, E_l + H_l + R_ll)})
eq_Rll = type('eq_Rll', (Equation,), {'__doc__': """R_ll as function of T_l and T_w.

    (Eq. 2 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(R_ll, a_sh*epsilon_l*sigm*(T_l**4 - T_w**4))})
eq_Hl = type('eq_Hl', (Equation,), {'__doc__': """H_l as function of T_l.

    (Eq. 3 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(H_l, a_sh*h_c*(-T_a + T_l))})
eq_El = type('eq_El', (Equation,), {'__doc__': """E_l as function of E_lmol.

    (Eq. 4 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(E_l, E_lmol*M_w*lambda_E)})
eq_Elmol = type('eq_Elmol', (Equation,), {'__doc__': """E_lmol as functino of g_tw and C_wl.

    (Eq. 5 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(E_lmol, g_tw*(-C_wa + C_wl))})
eq_gtw = type('eq_gtw', (Equation,), {'__doc__': """g_tw from g_sw and g_bw.

    (Eq. 6 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(g_tw, 1/(1/g_sw + 1/g_bw))})
eq_gbw_hc = type('eq_gbw_hc', (Equation,), {'__doc__': """g_bw as function of h_c.

    (Eq. B2 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(g_bw, a_s*h_c/(Le**(2/3)*c_pa*rho_a))})
eq_Cwl = type('eq_Cwl', (Equation,), {'__doc__': """C_wl as function of P_wl and T_l.

    (Eq. B4 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(C_wl, P_wl/(R_mol*T_l))})
eq_Pwl = type('eq_Pwl', (Equation,), {'__doc__': """Clausius-Clapeyron P_wl as function of T_l.

    (Eq. B3 in :cite:`hartmann_global_1994`)
    """, 'expr': Eq(P_wl, p_CC1*exp(-M_w*lambda_E*(-1/p_CC2 + 1/T_l)/R_mol))})
eq_Elmol_conv = type('eq_Elmol_conv', (Equation,), {'__doc__': """E_lmol as function of g_twmol and P_wl.

    (Eq. B6 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(E_lmol, g_twmol*(-P_wa + P_wl)/P_a)})
eq_gtwmol_gtw = type('eq_gtwmol_gtw', (Equation,), {'__doc__': """g_twmol as a function of g_tw.

    It uses eq_Elmol, eq_Cwl and eq_Elmol_conv.
    """, 'expr': Eq(g_twmol, g_tw*(P_a*P_wa*T_l - P_a*P_wl*T_a)/(R_mol*T_a*T_l*(P_wa - P_wl)))})
eq_gtwmol_gtw_iso = type('eq_gtwmol_gtw_iso', (Equation,), {'__doc__': """g_twmol as a function of g_tw at isothermal conditions.""", 'expr': Eq(g_twmol, P_a*g_tw/(R_mol*T_a))})
eq_hc = type('eq_hc', (Equation,), {'__doc__': """h_c as a function of Nu and L_l.

    (Eq. B10 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(h_c, Nu*k_a/L_l)})
eq_Re = type('eq_Re', (Equation,), {'__doc__': """Re as a function of v_w and L_l.

    (Eq. B11 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(Re, L_l*v_w/nu_a)})
eq_Gr = type('eq_Gr', (Equation,), {'__doc__': """Gr as function of air density within and outside of leaf.

    (Eq. B12 in :cite:`schymanski_leaf-scale_2017`)
    """, 'expr': Eq(Gr, L_l**3*g*(rho_a - rho_al)/(nu_a**2*rho_al))})
eq_Ew_fu = type('eq_Ew_fu', (Equation,), {'__doc__': """Penman water vapour diffusion.""", 'expr': Eq(E_w, f_u*(-P_wa + P_wl))})
eq_beta_B = type('eq_beta_B', (Equation,), {'__doc__': """Penman Bowen ratio.""", 'expr': Eq(beta_B, gamma_v*(-T_a + T_l)/(-P_wa + P_wl))})
eq_Penman_ass = type('eq_Penman_ass', (Equation,), {'__doc__': """""", 'expr': Eq(Delta_eTa, (-P_was + P_wl)/(-T_a + T_l))})
eq_Pwas_Ta = type('eq_Pwas_Ta', (Equation,), {'__doc__': """""", 'expr': Eq(P_was, p_CC1*exp(-M_w*lambda_E*(-1/p_CC2 + 1/T_a)/R_mol))})
eq_Deltaeta_T = type('eq_Deltaeta_T', (Equation,), {'__doc__': """""", 'expr': Eq(Delta_eTa, M_w*lambda_E*p_CC1*exp(-M_w*lambda_E*(-1/p_CC2 + 1/T_a)/R_mol)/(R_mol*T_a**2))})
eq_betaB_Pwas = type('eq_betaB_Pwas', (Equation,), {'__doc__': """""", 'expr': Eq(beta_B, gamma_v*(P_was - P_wl)/(Delta_eTa*(P_wa - P_wl)))})
eq_Ew_betaB = type('eq_Ew_betaB', (Equation,), {'__doc__': """""", 'expr': Eq(E_w, (-R_ll + R_s)/(beta_B + 1))})
eq_Pwl_fu = type('eq_Pwl_fu', (Equation,), {'__doc__': """""", 'expr': Eq(P_wl, (Delta_eTa*P_wa*f_u - Delta_eTa*R_ll + Delta_eTa*R_s + P_was*f_u*gamma_v)/(f_u*(Delta_eTa + gamma_v)))})
eq_Ew_P = type('eq_Ew_P', (Equation,), {'__doc__': """""", 'expr': Eq(E_w, (-Delta_eTa*R_ll + Delta_eTa*R_s - P_wa*f_u*gamma_v + P_was*f_u*gamma_v)/(Delta_eTa + gamma_v))})
eq_Tl_P = type('eq_Tl_P', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (Delta_eTa*T_a - P_was + P_wl)/Delta_eTa)})
eq_Tl_P = type('eq_Tl_P', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (Delta_eTa*T_a*f_u + P_wa*f_u - P_was*f_u - R_ll + R_s + T_a*f_u*gamma_v)/(f_u*(Delta_eTa + gamma_v)))})
eq_Tl_P1 = type('eq_Tl_P1', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, M_w*lambda_E*p_CC2/(M_w*lambda_E - R_mol*p_CC2*log(P_wl/p_CC1)))})
eq_El_fu_S = type('eq_El_fu_S', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, S*f_u*(-P_wa + P_wl))})
eq_Hl_Tl_P52 = type('eq_Hl_Tl_P52', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, f_u*gamma_v*(-T_a + T_l))})
eq_Hl_Pwl_P52 = type('eq_Hl_Pwl_P52', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, -f_u*gamma_v*(P_was - P_wl)/Delta_eTa)})
eq_El_P52 = type('eq_El_P52', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, S*(-Delta_eTa*(R_ll - R_s) - P_wa*f_u*gamma_v + P_was*f_u*gamma_v)/(Delta_eTa*S + gamma_v))})
eq_Hl_P52 = type('eq_Hl_P52', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, gamma_v*(P_wa*S*f_u - P_was*S*f_u - R_ll + R_s)/(Delta_eTa*S + gamma_v))})
eq_Pwl_P52 = type('eq_Pwl_P52', (Equation,), {'__doc__': """""", 'expr': Eq(P_wl, (Delta_eTa*P_wa*S*f_u - Delta_eTa*(R_ll - R_s) + P_was*f_u*gamma_v)/(f_u*(Delta_eTa*S + gamma_v)))})
eq_Tl_P52 = type('eq_Tl_P52', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (Delta_eTa*S*T_a*f_u + P_wa*S*f_u - P_was*S*f_u - R_ll + R_s + T_a*f_u*gamma_v)/(f_u*(Delta_eTa*S + gamma_v)))})
eq_S_gtwmol_fu = type('eq_S_gtwmol_fu', (Equation,), {'__doc__': """""", 'expr': Eq(S, M_w*g_twmol*lambda_E/(P_a*f_u))})
eq_gammav_hc_fu = type('eq_gammav_hc_fu', (Equation,), {'__doc__': """""", 'expr': Eq(gamma_v, a_sh*h_c/f_u)})
eq_Ew_conv = type('eq_Ew_conv', (Equation,), {'__doc__': """""", 'expr': Eq(E_w, -M_w*g_bw*lambda_E*(P_wa - P_wl)/(R_mol*T_a))})
eq_fu_gbw = type('eq_fu_gbw', (Equation,), {'__doc__': """""", 'expr': Eq(f_u, M_w*g_bw*lambda_E/(R_mol*T_a))})
eq_gammav_as = type('eq_gammav_as', (Equation,), {'__doc__': """""", 'expr': Eq(gamma_v, Le**0.666666666666667*R_mol*T_a*a_sh*c_pa*rho_a/(M_w*a_s*lambda_E))})
eq_S_gbw_gsw = type('eq_S_gbw_gsw', (Equation,), {'__doc__': """""", 'expr': Eq(S, g_sw/(g_bw + g_sw))})
eq_fu_ra_M = type('eq_fu_ra_M', (Equation,), {'__doc__': """""", 'expr': Eq(f_u, c_pa*rho_a/(gamma_v*r_a))})
eq_Ew_PM1 = type('eq_Ew_PM1', (Equation,), {'__doc__': """""", 'expr': Eq(E_w, (-Delta_eTa*R_ll + Delta_eTa*R_s - P_wa*c_pa*rho_a/r_a + P_was*c_pa*rho_a/r_a)/(Delta_eTa + gamma_v))})
eq_gammavs_M65 = type('eq_gammavs_M65', (Equation,), {'__doc__': """""", 'expr': Eq(gamma_v, gamma_v*(1 + r_s/r_a))})
eq_El_PM2 = type('eq_El_PM2', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, (-Delta_eTa*R_ll + Delta_eTa*R_s - P_wa*c_pa*rho_a/r_a + P_was*c_pa*rho_a/r_a)/(Delta_eTa + gamma_v*(1 + r_s/r_a)))})
eq_gammavs_MU = type('eq_gammavs_MU', (Equation,), {'__doc__': """""", 'expr': Eq(gamma_v, gamma_v*n_MU*(1 + r_s/r_a))})
eq_El_MU2 = type('eq_El_MU2', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, (-Delta_eTa*R_ll + Delta_eTa*R_s - P_wa*c_pa*rho_a/r_a + P_was*c_pa*rho_a/r_a)/(Delta_eTa + gamma_v*n_MU*(1 + r_s/r_a)))})
eq_gammav_MU = type('eq_gammav_MU', (Equation,), {'__doc__': """""", 'expr': Eq(gamma_v, P_a*c_pa/(epsilon*lambda_E))})
eq_epsilon = type('eq_epsilon', (Equation,), {'__doc__': """""", 'expr': Eq(epsilon, M_w*P_a/(R_mol*T_a*rho_a))})
eq_El_MU = type('eq_El_MU', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, a_s*epsilon*lambda_E*rho_a*(-P_wa + P_wl)/(P_a*(r_s + r_v)))})
eq_Hl_MU = type('eq_Hl_MU', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, a_sh*c_pa*rho_a*(-T_a + T_l)/r_a)})
eq_rv_gbw = type('eq_rv_gbw', (Equation,), {'__doc__': """""", 'expr': Eq(r_v, a_s/g_bw)})
eq_rs_gsw = type('eq_rs_gsw', (Equation,), {'__doc__': """""", 'expr': Eq(r_s, a_s/g_sw)})
eq_ra_hc = type('eq_ra_hc', (Equation,), {'__doc__': """""", 'expr': Eq(r_a, c_pa*rho_a/h_c)})
eq_El_cE = type('eq_El_cE', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, c_E*(-P_wa + P_wl))})
eq_Hl_cH = type('eq_Hl_cH', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, c_H*(-T_a + T_l))})
eq_gammav_cE = type('eq_gammav_cE', (Equation,), {'__doc__': """""", 'expr': Eq(gamma_v, c_H/c_E)})
eq_Hl_gammav = type('eq_Hl_gammav', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, c_E*gamma_v*(-T_a + T_l))})
eq_Hl_Pwl = type('eq_Hl_Pwl', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, -c_E*gamma_v*(P_was - P_wl)/Delta_eTa)})
eq_El_Delta = type('eq_El_Delta', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, (-Delta_eTa*(R_ll - R_s) - P_wa*c_H + P_was*c_H)/(Delta_eTa + c_H/c_E))})
eq_Hl_Delta = type('eq_Hl_Delta', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, c_H*(P_wa*c_E - P_was*c_E - R_ll + R_s)/(c_E*(Delta_eTa + c_H/c_E)))})
eq_Pwl_Delta = type('eq_Pwl_Delta', (Equation,), {'__doc__': """""", 'expr': Eq(P_wl, (Delta_eTa*P_wa*c_E - Delta_eTa*(R_ll - R_s) + P_was*c_H)/(c_E*(Delta_eTa + c_H/c_E)))})
eq_Tl_Delta = type('eq_Tl_Delta', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (-R_ll + R_s + T_a*c_H + c_E*(Delta_eTa*T_a + P_wa - P_was))/(Delta_eTa*c_E + c_H))})
eq_Tl_Delta1 = type('eq_Tl_Delta1', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, T_a + (-R_ll + R_s + c_E*(P_wa - P_was))/(Delta_eTa*c_E + c_H))})
eq_Twl_Delta2 = type('eq_Twl_Delta2', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, M_w*lambda_E/(R_mol*log(p_CC1*(Delta_eTa*c_E + c_H)*exp(M_w*lambda_E/(R_mol*p_CC2))/(Delta_eTa*P_wa*c_E - Delta_eTa*R_ll + Delta_eTa*R_s + P_was*c_H))))})
eq_ce_conv = type('eq_ce_conv', (Equation,), {'__doc__': """""", 'expr': Eq(c_E, M_w*g_twmol*lambda_E/P_a)})
eq_ch_hc = type('eq_ch_hc', (Equation,), {'__doc__': """""", 'expr': Eq(c_H, a_sh*h_c)})
eq_El_gammav = type('eq_El_gammav', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, -c_H*(P_wa - P_wl)/gamma_v)})
eq_El_Tl = type('eq_El_Tl', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, c_H*(-Delta_eTa*T_a + Delta_eTa*T_l - P_wa + P_was)/gamma_v)})
eq_El_Delta_a = type('eq_El_Delta_a', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, (-Delta_eTa*R_ll + Delta_eTa*R_s - P_wa*c_H + P_was*c_H)/(Delta_eTa + c_H/c_E))})
eq_Hl_Delta_a = type('eq_Hl_Delta_a', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, (P_wa*c_H - P_was*c_H - R_ll*c_H/c_E + R_s*c_H/c_E)/(Delta_eTa + c_H/c_E))})
eq_Tl_Delta_a = type('eq_Tl_Delta_a', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (Delta_eTa*T_a*c_H + P_wa*c_H - P_was*c_H + T_a*c_H**2/c_E - c_H*(R_ll - R_s)/c_E)/(c_H*(Delta_eTa + c_H/c_E)))})
eq_El_Delta_b = type('eq_El_Delta_b', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, -c_E*(Delta_eTa*R_ll - Delta_eTa*R_s + P_wa*c_H - P_was*c_H)/(Delta_eTa*c_E + c_H))})
eq_Hl_Delta_b = type('eq_Hl_Delta_b', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, c_H*(P_wa*c_E - P_was*c_E - R_ll + R_s)/(Delta_eTa*c_E + c_H))})
eq_Pwl_Delta_b = type('eq_Pwl_Delta_b', (Equation,), {'__doc__': """""", 'expr': Eq(P_wl, (Delta_eTa*P_wa*c_E - Delta_eTa*R_ll + Delta_eTa*R_s + P_was*c_H)/(Delta_eTa*c_E + c_H))})
eq_Tl_Delta_b = type('eq_Tl_Delta_b', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (P_wa*c_E - R_ll + R_s + T_a*c_H + c_E*(Delta_eTa*T_a - P_was))/(Delta_eTa*c_E + c_H))})
eq_El_Delta_MUcorr = type('eq_El_Delta_MUcorr', (Equation,), {'__doc__': """Corrected MU-equation""", 'expr': Eq(E_l, -a_s*epsilon*lambda_E*(Delta_eTa*r_a*(R_ll - R_s) + P_wa*a_sh*c_pa*rho_a - P_was*a_sh*c_pa*rho_a)/(Delta_eTa*a_s*epsilon*lambda_E*r_a + P_a*a_sh*c_pa*(r_a + r_s)))})
eq_El_MU_corr = type('eq_El_MU_corr', (Equation,), {'__doc__': """Corrected MU-equation, equivalent to Eq. 22 in \citet{schymanski_leaf-scale_2018}""", 'expr': Eq(E_l, (-Delta_eTa*R_ll + Delta_eTa*R_s - P_wa*a_sh*c_pa*rho_a/r_a + P_was*a_sh*c_pa*rho_a/r_a)/(Delta_eTa + a_sh*gamma_v*(1 + r_s/r_a)/a_s))})
eq_Rll_tang = type('eq_Rll_tang', (Equation,), {'__doc__': """Linearised R_ll.""", 'expr': Eq(R_ll, a_sh*epsilon_l*sigm*(-3*T_a**4 + 4*T_a**3*T_l - T_w**4))})
eq_El_Delta_Rlllin = type('eq_El_Delta_Rlllin', (Equation,), {'__doc__': """""", 'expr': Eq(E_l, c_E*(Delta_eTa*R_s - Delta_eTa*T_a**4*a_sh*epsilon_l*sigm + Delta_eTa*T_w**4*a_sh*epsilon_l*sigm - 4*P_wa*T_a**3*a_sh*epsilon_l*sigm - P_wa*c_H + 4*P_was*T_a**3*a_sh*epsilon_l*sigm + P_was*c_H)/(Delta_eTa*c_E + 4*T_a**3*a_sh*epsilon_l*sigm + c_H))})
eq_Hl_Delta_Rlllin = type('eq_Hl_Delta_Rlllin', (Equation,), {'__doc__': """""", 'expr': Eq(H_l, c_H*(P_wa*c_E - P_was*c_E + R_s - T_a**4*a_sh*epsilon_l*sigm + T_w**4*a_sh*epsilon_l*sigm)/(Delta_eTa*c_E + 4*T_a**3*a_sh*epsilon_l*sigm + c_H))})
eq_Pwl_Delta_Rlllin = type('eq_Pwl_Delta_Rlllin', (Equation,), {'__doc__': """""", 'expr': Eq(P_wl, (Delta_eTa*P_wa*c_E + Delta_eTa*R_s - Delta_eTa*T_a**4*a_sh*epsilon_l*sigm + Delta_eTa*T_w**4*a_sh*epsilon_l*sigm + 4*P_was*T_a**3*a_sh*epsilon_l*sigm + P_was*c_H)/(Delta_eTa*c_E + 4*T_a**3*a_sh*epsilon_l*sigm + c_H))})
eq_Tl_Delta_Rlllin = type('eq_Tl_Delta_Rlllin', (Equation,), {'__doc__': """""", 'expr': Eq(T_l, (P_wa*c_E + R_s + 3*T_a**4*a_sh*epsilon_l*sigm + T_a*c_H + T_w**4*a_sh*epsilon_l*sigm + c_E*(Delta_eTa*T_a - P_was))/(Delta_eTa*c_E + 4*T_a**3*a_sh*epsilon_l*sigm + c_H))})
