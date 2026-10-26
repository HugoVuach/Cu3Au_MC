"""Bragg-Williams à quatre sous-réseaux : approximation, pas résultat MC."""
import numpy as np
from scipy.optimize import brentq
from scipy.special import xlogy
from .model import KB

A = .0674625


def g(eta):
    p0, pf = (1+3*np.asarray(eta))/4, (1-np.asarray(eta))/4
    h = lambda p: xlogy(p, p)+xlogy(1-p, 1-p)
    return .25*h(p0)+.75*h(pf)


def gp(eta):
    return 3/16*np.log(((1+3*eta)*(3+eta))/(3*(1-eta)**2))


def coexistence():
    eta = brentq(lambda x: x*gp(x)-2*(g(x)-g(0)), .01, .999999)
    return {'eta_jump': float(eta), 'T_MF_K': float(2*A*eta/gp(eta)/KB),
            'T_spinodal_K': float(2*A/KB)}
