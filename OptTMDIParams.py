# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:38:17 2019

@author: Saptarshi
"""
from sympy import *
from sympy import symbols, diff, Matrix
mu, nu, b, w, wt, zd, phi = symbols('\mu, w_r, b, \omega, \omega_t, \zeta_d, \phi', positive=True, real=True)
init_printing()
M = Matrix([[1 + mu + b*(1-phi)**2, mu + b*(1-phi)], [mu + b*(1-phi), mu+b]])
K = Matrix([[wt**2, 0], [0, mu*nu**2*wt**2]])
C = Matrix([[0, 0], [0, 2*mu*nu*wt*zd]])

A = (-w**2*M + 1*I*w*C + K).inv()

G1 = A[0]

G11 = G1*conjugate(G1)
n,d = fraction(G11)
n1, d1 = fraction(G1)
"""
collected_numen = collect(expand(n), w)
collected_denom = collect(expand(d1), w)

b2 = simplify(collected_numen.coeff(w, 4))
b1 = simplify(collected_numen.coeff(w, 2))
b0 = simplify(collected_numen.coeff(w, 0))

a4 = simplify(collected_denom.coeff(w, 4))
a3 = simplify(collected_denom.coeff(w, 3))
a2 = simplify(collected_denom.coeff(w, 2))
a1 = simplify(collected_denom.coeff(w, 1))
a0 = simplify(collected_denom.coeff(w, 0))

d4 = a4
d3 = a3/(-1*I)
d2 = a2/(-1)
d1 = a1/(1*I)
d0 = a0

I4Num = Matrix([[0, b2, b1, b0], [-d4, d2, -d0, 0], [0, -d3, d1, 0], [0, d4, -d2, d0]])
I4Dem = Matrix([[d3, -d1, 0, 0], [-d4, d2, -d0, 0], [0, -d3, d1, 0], [0, d4, -d2, d0]])

Integration = (pi*I4Num.det())/(d4*I4Dem.det())

Int  = simplify(Integration)

# Divison here done manually
OptTuning = simplify(diff(Int, nu))
OptTun  = OptTuning/(pi/fraction(OptTuning)[1])

OptTun  = collect(OptTun, nu)
OptTunC = simplify(OptTun.coeff(nu, 0))
OptTunB = simplify(OptTun.coeff(nu, 2))
OptTunA = simplify(OptTun.coeff(nu, 4))

# Divison here done manually
OptDamp = simplify(diff(Int, zd))
OptDmp  = OptDamp/(pi/fraction(OptDamp)[1])
OptDmp  = collect(OptDmp, zd)
OptDmpC = simplify(OptDmp.coeff(zd, 0))
OptDmpB = simplify(OptDmp.coeff(zd, 1))
OptDmpA = simplify(OptDmp.coeff(zd, 2))
soloptdamp = solve(OptDmp, zd)
"""