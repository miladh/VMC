# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:08:08 2013

@author: Milad H. Mobarhan
"""
from sympy import (Symbol,
                   symbols,
                   sqrt,
                   exp,
                   diff,
                   printing,
                   sympify)

import numpy as np

r = Symbol('r', real=True, positive=True)
R = Symbol('R', real=True, positive=True)
rp1Sym, rp2Sym, rp1Sym_2, rp2Sym_2 = symbols('rp1 rp2 rp1^2 rp2^2', real=True, positive=True)
xp1, xp2, yp1, yp2, zp1, zp2 = symbols('xp1 xp2 yp1 yp2 zp1 zp2', real=True)


a = Symbol('a', real=True, positive=True)
x, y, z = symbols('x y z', real=True)
rNorm = sqrt(x**2 + y**2 + z**2)

X, Y, Z = symbols('X Y Z', real=True)

RNorm = sqrt(X**2 + Y**2 + Z**2)

rp1_2 = (x + X)**2 + (y+Y)**2 + (z+Z)**2
rp1 = sqrt(rp1_2)

rp2_2 = (x - X)**2 + (y-Y)**2 + (z-Z)**2
rp2 = sqrt(rp2_2)

def expFactor(rp):
    return exp(-a*rp)

phi = expFactor(rp1) + expFactor(rp2)

xyz = [x, y, z]

grad = []
for xi in xyz:
    grad.append(diff(phi, xi))

lapl = 0
for xi, dell in zip(xyz, grad):
    print xi, ": ", dell.subs(rp1, rp1Sym).subs(rp2, rp2Sym).subs(X+ x, xp1).subs(-X + x, xp2).\
                                                            subs(Y+ y, yp1).subs(-Y + y, yp2).\
                                                            subs(Z+ z, zp1).subs(-Z + z, zp2)
    lapl += diff(dell, xi)
#lapl = lapl.collect(a)
lapl = lapl.subs(rp1, rp1Sym).subs(rp2, rp2Sym)

p2Terms = []
p1Terms = []

for term in lapl.as_ordered_terms():
    if exp(-a*rp1Sym) in term:
        p1Terms.append(term)
    elif exp(-a*rp2Sym) in term:
        p2Terms.append(term)
    else:
        print "you fucked up...", term

def getLapl(rp, terms):
    lapl = sum(terms)
    lapl = lapl.collect(expFactor(rp))/expFactor(rp)
    #for term1, term2 in zip(p1Terms, p2Terms):
    #    print term1, "   |   " ,term2
    print
    lapl = lapl.collect(a/rp).subs(rp2_2, rp2Sym**2).subs(rp1_2, rp1Sym**2)
    lapl = sympify(str(lapl).replace("-(-X - x)*(X + x) - (-Y - y)*(Y + y) - (-Z - z)*(Z + z)", "rp1**2"))
    lapl = sympify(str(lapl).replace("(X - x)**2 - (-Y + y)*(Y - y) - (-Z + z)*(Z - z)", "rp2**2"))
    #lapl1 = lapl1.factor(a)
    return lapl*expFactor(rp)
    
lapl = getLapl(rp1Sym, p1Terms) + getLapl(rp2Sym, p2Terms)

print lapl
