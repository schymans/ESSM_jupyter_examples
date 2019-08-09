from essm.variables._core import BaseVariable, Variable
from essm.equations import Equation
from sympy import Abs, Derivative, Eq, exp, Integral, log, Piecewise, sqrt
from sympy.physics.units import joule, kilogram, mole, pascal, kelvin, meter, watt, second
c_pi = type('c_pi', (Variable,), {'__doc__': """Heat capacity of insulation material.""", 'unit': joule/(kelvin*kilogram), 'assumptions': {'real': True},         'latex_name': r'c_{pi}', 'default': None, 'expr': None})
lambda_i = type('lambda_i', (Variable,), {'__doc__': """Heat conductivity of insulation material.""", 'unit': joule/(kelvin*meter*second), 'assumptions': {'real': True},         'latex_name': r'lambda_i', 'default': None, 'expr': None})
rho_i = type('rho_i', (Variable,), {'__doc__': """Density of insulation material.""", 'unit': kilogram/meter**3, 'assumptions': {'real': True},         'latex_name': r'rho_i', 'default': None, 'expr': None})
L_i = type('L_i', (Variable,), {'__doc__': """Thickness of insulation material.""", 'unit': meter, 'assumptions': {'real': True},         'latex_name': r'L_i', 'default': None, 'expr': None})
A_i = type('A_i', (Variable,), {'__doc__': """Conducting area of insulation material.""", 'unit': meter**2, 'assumptions': {'real': True},         'latex_name': r'A_i', 'default': None, 'expr': None})
Q_i = type('Q_i', (Variable,), {'__doc__': """Heat conduction through insulation material.""", 'unit': joule/second, 'assumptions': {'real': True},         'latex_name': r'Q_i', 'default': None, 'expr': None})
dT_i = type('dT_i', (Variable,), {'__doc__': """Temperature increment of insulation material.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'dT_i', 'default': None, 'expr': None})
alpha_a = type('alpha_a', (Variable,), {'__doc__': """Thermal diffusivity of dry air.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'\alpha_a', 'default': None, 'expr': None})
c_pa = type('c_pa', (Variable,), {'__doc__': """Specific heat of dry air.""", 'unit': joule/(kelvin*kilogram), 'assumptions': {'real': True},         'latex_name': r'c_{pa}', 'default': 1010.0, 'expr': None})
c_pamol = type('c_pamol', (Variable,), {'__doc__': """Molar specific heat of dry air.      https://en.wikipedia.org/wiki/Heat_capacity#Specific_heat_capacity     """, 'unit': joule/(kelvin*mole), 'assumptions': {'real': True},         'latex_name': r'c_{pa,mol}', 'default': 29.19, 'expr': None})
c_pv = type('c_pv', (Variable,), {'__doc__': """Specific heat of water vapour at 300 K.      http://www.engineeringtoolbox.com/water-vapor-d_979.html     """, 'unit': joule/(kelvin*kilogram), 'assumptions': {'real': True},         'latex_name': r'c_{pv}', 'default': 1864, 'expr': None})
C_wa = type('C_wa', (Variable,), {'__doc__': """Concentration of water in air.""", 'unit': mole/meter**3, 'assumptions': {'real': True},         'latex_name': r'C_{wa}', 'default': None, 'expr': None})
D_va = type('D_va', (Variable,), {'__doc__': """Binary diffusion coefficient of water vapour in air.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'D_{va}', 'default': None, 'expr': None})
g = type('g', (Variable,), {'__doc__': """Gravitational acceleration.""", 'unit': meter/second**2, 'assumptions': {'real': True},         'latex_name': r'g', 'default': 9.81, 'expr': None})
Gr = type('Gr', (Variable,), {'__doc__': """Grashof number.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Gr_L}', 'default': None, 'expr': None})
h_c = type('h_c', (Variable,), {'__doc__': """Average 1-sided convective heat transfer coefficient.""", 'unit': joule/(kelvin*meter**2*second), 'assumptions': {'real': True},         'latex_name': r'h_c', 'default': None, 'expr': None})
k_a = type('k_a', (Variable,), {'__doc__': """Thermal conductivity of dry air.""", 'unit': joule/(kelvin*meter*second), 'assumptions': {'real': True},         'latex_name': r'k_a', 'default': None, 'expr': None})
lambda_E = type('lambda_E', (Variable,), {'__doc__': """Latent heat of evaporation.""", 'unit': joule/kilogram, 'assumptions': {'real': True},         'latex_name': r'\lambda_E', 'default': 2450000.0, 'expr': None})
Le = type('Le', (Variable,), {'__doc__': """Lewis number.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Le}', 'default': None, 'expr': None})
M_air = type('M_air', (Variable,), {'__doc__': """Molar mass of air.      http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html     """, 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_{air}', 'default': 0.02897, 'expr': None})
M_N2 = type('M_N2', (Variable,), {'__doc__': """Molar mass of nitrogen.""", 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_{N_2}', 'default': 0.028, 'expr': None})
M_O2 = type('M_O2', (Variable,), {'__doc__': """Molar mass of oxygen.""", 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_{O_2}', 'default': 0.032, 'expr': None})
M_w = type('M_w', (Variable,), {'__doc__': """Molar mass of water.""", 'unit': kilogram/mole, 'assumptions': {'real': True},         'latex_name': r'M_w', 'default': 0.018, 'expr': None})
nu_a = type('nu_a', (Variable,), {'__doc__': """Kinematic viscosity of dry air.""", 'unit': meter**2/second, 'assumptions': {'real': True},         'latex_name': r'\nu_a', 'default': None, 'expr': None})
Nu = type('Nu', (Variable,), {'__doc__': """Average Nusselt number over given length.""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Nu_L}', 'default': None, 'expr': None})
P_a = type('P_a', (Variable,), {'__doc__': """Air pressure.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_a', 'default': None, 'expr': None})
Pr = type('Pr', (Variable,), {'__doc__': """Prandtl number (0.71 for air).""", 'unit': 1, 'assumptions': {'real': True},         'latex_name': r'N_{Pr}', 'default': None, 'expr': None})
P_N2 = type('P_N2', (Variable,), {'__doc__': """Partial pressure of nitrogen.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{N2}', 'default': None, 'expr': None})
P_O2 = type('P_O2', (Variable,), {'__doc__': """Partial pressure of oxygen.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{O2}', 'default': None, 'expr': None})
P_wa = type('P_wa', (Variable,), {'__doc__': """Partial pressure of water vapour in air.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{wa}', 'default': None, 'expr': None})
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
P_g = type('P_g', (Variable,), {'__doc__': """Pressure of gas.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_g', 'default': None, 'expr': None})
V_g = type('V_g', (Variable,), {'__doc__': """Volume of gas.""", 'unit': meter**3, 'assumptions': {'real': True},         'latex_name': r'V_g', 'default': None, 'expr': None})
n_g = type('n_g', (Variable,), {'__doc__': """Amount of gas.""", 'unit': mole, 'assumptions': {'real': True},         'latex_name': r'n_g', 'default': None, 'expr': None})
n_w = type('n_w', (Variable,), {'__doc__': """Amount of water.""", 'unit': mole, 'assumptions': {'real': True},         'latex_name': r'n_w', 'default': None, 'expr': None})
T_g = type('T_g', (Variable,), {'__doc__': """Temperature of gas.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_g', 'default': None, 'expr': None})
Delta_Pwa = type('Delta_Pwa', (Variable,), {'__doc__': """Slope of saturated vapour pressure, $\partial P_{wa} / \partial T_g$""", 'unit': pascal/kelvin, 'assumptions': {'real': True},         'latex_name': r'\Delta', 'default': None, 'expr': Derivative(P_wa, T_g)})
x = type('x', (Variable,), {'__doc__': """Positive real variable.""", 'unit': 1, 'assumptions': {'positive': True, 'real': True},         'latex_name': r'x', 'default': None, 'expr': None})
p_CC1 = type('p_CC1', (Variable,), {'__doc__': """Internal parameter of eq_Pwl.""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'611', 'default': 611.0, 'expr': None})
p_CC2 = type('p_CC2', (Variable,), {'__doc__': """Internal parameter of eq_Pwl.""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'273', 'default': 273.0, 'expr': None})
T_a1 = type('T_a1', (Variable,), {'__doc__': """Air temperature""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_{a1}', 'default': None, 'expr': None})
T_a2 = type('T_a2', (Variable,), {'__doc__': """Air temperature""", 'unit': kelvin, 'assumptions': {'real': True},         'latex_name': r'T_{a2}', 'default': None, 'expr': None})
P_wa1 = type('P_wa1', (Variable,), {'__doc__': """P_wa at T1""", 'unit': pascal, 'assumptions': {'real': True},         'latex_name': r'P_{wa1}', 'default': None, 'expr': None})
eq_Qi = type('eq_Qi', (Equation,), {'__doc__': """Calculate ....      :cite:`schymanski_leaf-scale_2017`     """, 'expr': Eq(Q_i, A_i*dT_i*lambda_i/L_i)})
eq_Le = type('eq_Le', (Equation,), {'__doc__': """Le as function of alpha_a and D_va.      (Eq. B3 in :cite:`schymanski_leaf-scale_2017`)     """, 'expr': Eq(Le, alpha_a/D_va)})
eq_Cwa = type('eq_Cwa', (Equation,), {'__doc__': """C_wa as a function of P_wa and T_a.      (Eq. B9 in :cite:`schymanski_leaf-scale_2017`)     """, 'expr': Eq(C_wa, P_wa/(R_mol*T_a))})
eq_Nu_forced_all = type('eq_Nu_forced_all', (Equation,), {'__doc__': """Nu as function of Re and Re_c under forced conditions.      (Eqs. B13--B15 in :cite:`schymanski_leaf-scale_2017`)     """, 'expr': Eq(Nu, -Pr**(1/3)*(-37*Re**(4/5) + 37*(Re + Re_c - Abs(Re - Re_c)/2)**(4/5) - 664*sqrt(Re + Re_c - Abs(Re - Re_c)/2))/1000)})
eq_Dva = type('eq_Dva', (Equation,), {'__doc__': """D_va as a function of air temperature.      (Table A.3 in :cite:`monteith_principles_2007`)     """, 'expr': Eq(D_va, T_a*p_Dva1 - p_Dva2)})
eq_alphaa = type('eq_alphaa', (Equation,), {'__doc__': """alpha_a as a function of air temperature.      (Table A.3 in :cite:`monteith_principles_2007`)     """, 'expr': Eq(alpha_a, T_a*p_alpha1 - p_alpha2)})
eq_ka = type('eq_ka', (Equation,), {'__doc__': """k_a as a function of air temperature.      (Table A.3 in :cite:`monteith_principles_2007`)     """, 'expr': Eq(k_a, T_a*p_ka1 + p_ka2)})
eq_nua = type('eq_nua', (Equation,), {'__doc__': """nu_a as a function of air temperature.      (Table A.3 in :cite:`monteith_principles_2007`)     """, 'expr': Eq(nu_a, T_a*p_nua1 - p_nua2)})
eq_rhoa_Pwa_Ta = type('eq_rhoa_Pwa_Ta', (Equation,), {'__doc__': """rho_a as a function of P_wa and T_a.      (Eq. B20 in :cite:`schymanski_leaf-scale_2017`)     """, 'expr': Eq(rho_a, (M_N2*P_N2 + M_O2*P_O2 + M_w*P_wa)/(R_mol*T_a))})
eq_Pa = type('eq_Pa', (Equation,), {'__doc__': """Calculate air pressure.      From partial pressures of N2, O2 and H2O, following Dalton's law of     partial pressures.     """, 'expr': Eq(P_a, P_N2 + P_O2 + P_wa)})
eq_PN2_PO2 = type('eq_PN2_PO2', (Equation,), {'__doc__': """Calculate P_N2 as a function of P_O2.      It follows Dalton's law of partial pressures.     """, 'expr': Eq(P_N2, P_O2*x_N2/x_O2)})
eq_PO2 = type('eq_PO2', (Equation,), {'__doc__': """Calculate P_O2 as a function of P_a, P_N2 and P_wa.""", 'expr': Eq(P_O2, (P_a*x_O2 - P_wa*x_O2)/(x_N2 + x_O2))})
eq_PN2 = type('eq_PN2', (Equation,), {'__doc__': """Calculate P_N2 as a function of P_a, P_O2 and P_wa.""", 'expr': Eq(P_N2, (P_a*x_N2 - P_wa*x_N2)/(x_N2 + x_O2))})
eq_rhoa = type('eq_rhoa', (Equation,), {'__doc__': """Calculate rho_a from T_a, P_a and P_wa.""", 'expr': Eq(rho_a, (x_N2*(M_N2*P_a - P_wa*(M_N2 - M_w)) + x_O2*(M_O2*P_a - P_wa*(M_O2 - M_w)))/(R_mol*T_a*x_N2 + R_mol*T_a*x_O2))})
eq_ideal_gas_law = type('eq_ideal_gas_law', (Equation,), {'__doc__': """Ideal gas law.""", 'expr': Eq(P_g*V_g, R_mol*T_g*n_g)})
eq_Pg = type('eq_Pg', (Equation,), {'__doc__': """Calculate pressure of ideal gas.""", 'expr': Eq(P_g, R_mol*T_g*n_g/V_g)})
eq_Pwa_nw = type('eq_Pwa_nw', (Equation,), {'__doc__': """Calculate vapour pressure from amount of water in gas.""", 'expr': Eq(P_wa, R_mol*T_g*n_w/V_g)})
eq_Pwa_CC = type('eq_Pwa_CC', (Equation,), {'__doc__': """Clausius-Clapeyron P_wa as function of T_g.       Eq. B3 in :cite{hartmann_global_1994}     """, 'expr': Eq(P_wa, p_CC1*exp(-M_w*lambda_E*(-1/p_CC2 + 1/T_g)/R_mol))})
eq1 = type('eq1', (Equation,), {'__doc__': """Test""", 'expr': Eq(P_wa, Piecewise((0, T_a < 0), (p_CC1*exp(-M_w*lambda_E*(-1/p_CC2 + 1/T_g)/R_mol), True)))})
eq_Pwa_Delta = type('eq_Pwa_Delta', (Equation,), {'__doc__': """P_wa deduced from the integral of Delta""", 'expr': Eq(P_wa, P_wa1 + Integral(Delta_Pwa, (T_g, T_a1, T_a2)))})