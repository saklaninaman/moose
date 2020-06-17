#!/usr/bin/env python3
#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def expected(x):
    ee = 1.2
    nu = 0.3
    t = 0.0002
    ll = 10
    z = 0.5
    c = 0.5
    gg = ee / 2 / (1+nu)
    beta = 3 * t * (1 - nu * nu) / 4 / c / c / c / ee
    dd = - beta * nu / (1 - nu)
    delta = beta / 3
    gamma = -3 * t / 4 / c / gg
    alpha = dd / 3 + t / 4 / c / c / c / gg
    ux = beta * x * z * (2 * ll - x) + alpha * pow(z, 3)
    uz = delta * x * x * (x - 3 * ll) + gamma * x + dd * z * z * (ll - x)
    return (ux, uz)

def expected2(x):
    ee = 1.2
    nu = 0.3
    aa = 1.11E-2
    ll = 10
    z = 0.5
    c = 0.5
    y = 0
    gg = ee / 2 / (1+nu)
    dd = -nu * aa
    ux = aa * x * z
    uy = dd * z * y
    uz = -0.5 * aa * x * x + 0.5 * dd * (z * z - y * y)
    return (ux, uy, uz)

def expected_slippery(x):
    ee = 1.2
    nu = 0.3
    ll = 10
    h = 0.5
    s = -2E-4
    gg = ee / 2 / (1 + nu)
    bb = ee * h * h / 12 / (1 - nu * nu)
    thy = 0.5 * s * x * (x - 2 * ll) / bb
    uz = 2 * s * x / gg + 2 * s * (1 - nu * nu) * x * x * (3 * ll - x) / ee / h / h
    return (thy, uz)

def expected_slippery_h(h):
    ee = 1.2
    nu = 0.3
    ll = 10
    x = 10
    s = -2E-4
    gg = ee / 2 / (1 + nu)
    bb = ee * h * h / 12 / (1 - nu * nu)
    thy = 0.5 * s * x * (x - 2 * ll) / bb
    uz = 2 * s * x / gg + 2 * s * (1 - nu * nu) * x * x * (3 * ll - x) / ee / h / h
    return (thy, uz)

def solid_bar():
    f = open("../../tests/static_deformations/gold/beam_cosserat_01_soln_0001.csv")
    data = [map(float, line.strip().split(",")) for line in f.readlines()[1:12]]
    f.close()
    return data

def solid_bar2():
    f = open("../../tests/static_deformations/gold/beam_cosserat_02_apply_stress_soln_0001.csv")
    data = [map(float, line.strip().split(",")) for line in f.readlines()[1:12]]
    f.close()
    return data

def solid_bar_slippery():
    f = open("../../tests/static_deformations/gold/beam_cosserat_01_slippery_soln_0001.csv")
    data = [map(float, line.strip().split(",")) for line in f.readlines()[1:12]]
    f.close()
    return data

def solid_bar_slippery_h():
    # these data were generated by hand using beam_cosserat_01_slippery.i with different values of layer_thickness (and nx=800)
    data = [(0.1, -60.3), (0.2, -15.1), (0.3, -6.74), (0.4, -3.8), (0.5, -2.4), (0.9, -0.76)]
    return data



xpoints = np.arange(0, 10.05, 0.1)
hpoints = np.arange(0.09, 1, 0.01)
moosex = [i for i in range(11)]
moose = solid_bar()
moose2 = solid_bar2()
moose_slippery = solid_bar_slippery()
mooseh = [0.1, 0.2, 0.3, 0.4, 0.5, 0.9]
moose_slippery_h = solid_bar_slippery_h()

plt.figure()
plt.plot(xpoints, expected(xpoints)[0], 'k-', linewidth = 1.0, label = 'expected u_x')
plt.plot(xpoints, expected(xpoints)[1], 'r-', linewidth = 1.0, label = 'expected u_z')
plt.plot(moosex, [d[4] for d in moose], 'ks', markersize = 10.0, label = 'MOOSE disp_x')
plt.plot(moosex, [d[5] for d in moose], 'r^', markersize = 10.0, label = 'MOOSE disp_z')
plt.legend(loc = 'lower left')
plt.xlabel("x (m)")
plt.ylabel("displacement (m)")
plt.title("Beam deformation")
#plt.savefig("cosserat_beam_disp.pdf")

plt.figure()
plt.plot(xpoints, expected2(xpoints)[0], 'k-', linewidth = 1.0, label = 'expected u_x')
plt.plot(xpoints, expected2(xpoints)[2], 'r-', linewidth = 1.0, label = 'expected u_z')
plt.plot(moosex, [d[9] for d in moose2], 'ks', markersize = 10.0, label = 'MOOSE disp_x')
plt.plot(moosex, [d[11] for d in moose2], 'r^', markersize = 10.0, label = 'MOOSE disp_z')
plt.legend(loc = 'lower left')
plt.xlabel("x (m)")
plt.ylabel("displacement (m)")
plt.title("Beam deformation")
#plt.savefig("cosserat_beam_disp_2.pdf")

plt.figure()
plt.plot(xpoints, expected_slippery(xpoints)[0], 'k-', linewidth = 1.0, label = 'expected th_y')
plt.plot(xpoints, expected_slippery(xpoints)[1], 'r-', linewidth = 1.0, label = 'expected u_z')
plt.plot(moosex, [d[11] for d in moose_slippery], 'ks', markersize = 10.0, label = 'MOOSE wc_y')
plt.plot(moosex, [d[5] for d in moose_slippery], 'r^', markersize = 10.0, label = 'MOOSE disp_z')
plt.legend(loc = 'lower left')
plt.xlabel("x (m)")
plt.ylabel("disp (m) and rot")
plt.title("Slippery beam deformation")
#plt.savefig("cosserat_beam_disp_slippery.pdf")

plt.figure()
plt.plot(hpoints, expected_slippery_h(hpoints)[1], 'k-', linewidth = 1.0, label = 'expected')
plt.plot(mooseh, [d[1] for d in moose_slippery_h], 'ks', markersize = 10.0, label = 'MOOSE')
plt.legend(loc = 'lower right')
plt.xlabel("h (m)")
plt.ylabel("deflection (m)")
plt.title("End-point deflection in slippery Cosserat bar")
plt.savefig("cosserat_beam_disp_slippery_h.pdf")

sys.exit(0)
