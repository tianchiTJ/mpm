import math

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.cm as cm

import json

# General paramters for plot
# Psi to Pa
psi_pa = 6894.75729
# Figure size
plt.rcParams['figure.figsize'] = (11.0, 9.0)
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['figure.dpi'] = 200
# label font
font_label = {'family': 'Times New Roman',
              'weight': '700',
              'style': 'italic',
              'size': 20,
              }
# legend
font_legend = {'family': 'Times New Roman',
               'weight': 'normal',
               'style': 'italic',
               'size': 8,
               }
# Colors
colors = cm.rainbow(np.linspace(0, 1, 5))
markers = ['o', 's', '^', 'x', '+']
c = colors[0]
marker = markers[0]
name = 'pconf = 500Psi'
######################################################
# Test data
# Number of test cases
triaxial_data = 15
# Path_16_x
path_16 = ["" for i in range(triaxial_data)]
for i in range(triaxial_data):
    path_16[i] = "./Exxon_data/16-"+str(i+1)+".txt"
# Number of simulation cases
simulation_data = 5
# path_simulation_x
path_simulation_p_q = ["" for i in range(simulation_data)]
path_simulation_strain = ["" for i in range(simulation_data)]
path_simulation_e = ["" for i in range(simulation_data)]
path_simulation_pc = ["" for i in range(simulation_data)]
path_simulation_r = ["" for i in range(simulation_data)]
for i in range(simulation_data):
    path_simulation_p_q[i] = "./Exxon_data/BMCC_subloading_results/MCC_bonded_subloading_drained_p_q_" + \
        str(i+1)+".txt"
    path_simulation_strain[i] = "./Exxon_data/BMCC_subloading_results/MCC_bonded_subloading_drained_strain_" + \
        str(i+1)+".txt"
    path_simulation_e[i] = "./Exxon_data/BMCC_subloading_results/MCC_bonded_subloading_drained_void_ratio_" + \
        str(i+1)+".txt"
    path_simulation_pc[i] = "./Exxon_data/BMCC_subloading_results/MCC_bonded_subloading_drained_pc_pcc_pcd_" + \
        str(i+1)+".txt"
    path_simulation_r[i] = "./Exxon_data/BMCC_subloading_results/MCC_bonded_subloading_drained_r_ratio_" + \
        str(i+1)+".txt"
######################################################
# Yield surface and critical lines
with open('./Exxon_data/BMCC-parameters.json') as json_file:
    data = json.load(json_file)
# pc0
pc0 = data['pc0']
# Critical slope
m = data['m']
# Bonded pcc and pcd
pcd0 = data['mc_a']*pow(data['s_h'], data['mc_b'])
pcc0 = data['mc_c']*pow(data['s_h'], data['mc_d'])
# Confining pressure
pconf = [500, 1000, 2500, 5000, 8000]
q_max = [4E7, 5E7, 8E7, 9E7, 1.1E8]
# Compute current pcc
# pcc_current = [0]*simulation_data
# for i in range(simulation_data):
#     pcc_current[i] = -(q_max[i] / 3 - q_max[i] / m + pconf[i]*psi_pa)
# Number of points
n = 1000000
# Critical line
q_critical = [0]*n
q_critical_bonded = [0]*n
mul = 1.1
# Yield surface
p_yield = [0]*n
q_yield = [0]*n
p_yield_bonded = [0]*n
q_yield_bonded = [0]*n
# Plot yield surface and critical lines
for i in range(n):
    p_yield[i] = pc0 / n * i
    q_yield[i] = math.sqrt(-m*m*(p_yield[i])*(p_yield[i]-pc0))
    q_critical[i] = m*p_yield[i]*mul
######################################################
# Figure 1 q-p curves
# plot
plt.figure()
figsize = 11, 9
figure, ax = plt.subplots(figsize=figsize)
# Axial ranges
# plt.axis([-np.amax(pcc_current), (pc0+np.amax(pcc_current)+pcd0)
#           * mul, 0, (pc0+np.amax(pcc_current)+pcd0)*mul*m])
plt.axis([-pcc0, (pc0+pcc0+pcd0)
          * mul, 0, (pc0+pcc0+pcd0)*mul*m])
plt.tick_params(labelsize=16)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
# xmajorFormatter = FormatStrFormatter('%1.1f')
# ax.xaxis.set_major_formatter(xmajorFormatter)

# pc pcc pcd
for i in range(simulation_data):
    p_simulation = np.loadtxt(
        path_simulation_p_q[i], delimiter='\t', usecols=(0))
    q_simulation = np.loadtxt(
        path_simulation_p_q[i], delimiter='\t', usecols=(1))
    pc_simulation = np.loadtxt(
        path_simulation_pc[i], delimiter='\t', usecols=(0))
    pcc_simulation = np.loadtxt(
        path_simulation_pc[i], delimiter='\t', usecols=(1))
    pcd_simulation = np.loadtxt(
        path_simulation_pc[i], delimiter='\t', usecols=(2))
    r_simulation = np.loadtxt(
        path_simulation_r[i], delimiter='\t', usecols=(0))
    pc0 = pc_simulation[0]
    pcc0 = pcc_simulation[0]
    pcd0 = pcd_simulation[0]
    r0 = r_simulation[0]
    # Set color
    if (p_simulation[0] > (7000*psi_pa)):
        c = colors[4]
        marker = markers[4]
        name = 'pconf = 8000Psi'
    elif (p_simulation[0] > (4000*psi_pa)):
        c = colors[3]
        marker = markers[3]
        name = 'pconf = 5000Psi'
    elif (p_simulation[0] > (2300*psi_pa)):
        c = colors[2]
        marker = markers[2]
        name = 'pconf = 2500Psi'
    elif (p_simulation[0] > (800*psi_pa)):
        c = colors[1]
        marker = markers[1]
        name = 'pconf = 1000Psi'
    else:
        c = colors[0]
        marker = markers[0]
        name = 'pconf = 500Psi'
    for j in range(len(pc_simulation)):
        if r_simulation[j]*(pc_simulation[j]+pcc_simulation[j]+pcd_simulation[j]) > r0*(pc0+pcc0+pcd0):
            pc0 = pc_simulation[j]
            pcc0 = pcc_simulation[j]
            pcd0 = pcd_simulation[j]
            r0 = r_simulation[j]
    # Plot bonded yield surface
    for k in range(n):
        p_yield_bonded[k] = r0*(pc0+pcc0*2+pcd0) / n * k - pcc0
        q_yield_bonded[k] = math.sqrt(-m*m*(p_yield_bonded[k]+pcc0)
                                      * (p_yield_bonded[k]-r0*(pc0+pcd0+pcc0)))
        q_critical_bonded[k] = m*p_yield_bonded[k]*mul + pcc0

    plt.plot(p_yield_bonded, q_yield_bonded, color=c, linestyle=':',
             label='Bonded yield surface'+str(i+1)+", "+name)

    # Plot simulation data
    plt.plot(p_simulation, q_simulation, color=c,
             alpha=1, label='BMCC Model'+str(i+1))

# plot bonded critical lines
# plt.plot([i * mul for i in p_yield_bonded],
#         q_critical_bonded, 'r-', alpha=1, label='bonded critical line')

# Plot yield surface
plt.plot(p_yield, q_yield, linewidth=3, color='k', label='Yield surface')
# Bonded yield surface
for i in range(n):
    p_yield_bonded[i] = (pc0+pcc0*2+pcd0) / n * i - pcc0
    q_yield_bonded[i] = math.sqrt(-m*m*(p_yield_bonded[i]+pcc0)
                                  * (p_yield_bonded[i]-pc0-pcd0-pcc0))
# -------------------------------------
plt.plot([i * mul for i in p_yield], q_critical,
         linewidth=3, color='k', label='Critical line')
# plot p-q
# Plot traxial test data
for i in range(triaxial_data):
    # p_16_x
    p_16 = psi_pa * np.loadtxt(path_16[i], delimiter='\t',
                               skiprows=2, usecols=(0))
    # q_16_x
    q_16 = psi_pa * np.loadtxt(path_16[i], delimiter='\t',
                               skiprows=2, usecols=(1))
    # Set color
    if (p_16[0] > (7000*psi_pa)):
        c = colors[4]
        marker = markers[4]
        name = 'pconf = 8000Psi'
    elif (p_16[0] > (4000*psi_pa)):
        c = colors[3]
        marker = markers[3]
        name = 'pconf = 5000Psi'
    elif (p_16[0] > (2300*psi_pa)):
        c = colors[2]
        marker = markers[2]
        name = 'pconf = 2500Psi'
    elif (p_16[0] > (800*psi_pa)):
        c = colors[1]
        marker = markers[1]
        name = 'pconf = 1000Psi'
    else:
        c = colors[0]
        marker = markers[0]
        name = 'pconf = 500Psi'
    # plot data
    plt.scatter(p_16, q_16, color=c, s=20, marker=marker,
                alpha=1, label='Traxial test data 16-'+str(i+1)+", "+name)

legend = plt.legend(prop=font_legend, loc='upper left')
# axis label name
plt.xlabel('p\' / Pa', font_label)
plt.ylabel('q / Pa', font_label)
# Save figure
plt.savefig('p-q.png')
######################################################
# Figure 2 q-astrain curves
# Plot
plt.figure()
figsize = 11, 9
figure, ax = plt.subplots(figsize=figsize)
# Axial accuracy
x_ticks = np.arange(0, 0.06, 0.01)
y_ticks = np.arange(0, 120000000, 10000000)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
# Axial ranges
plt.axis([0, 0.06, 0, 120000000])
plt.tick_params(labelsize=16)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
# axis label name
plt.xlabel('axial strain', font_label)
plt.ylabel('q / Pa', font_label)
# Plot traxial test data
for i in range(triaxial_data):
    # p_16_x
    p_16 = psi_pa * np.loadtxt(path_16[i], delimiter='\t',
                               skiprows=2, usecols=(0))
    # astrain_16_x
    astrain_16 = np.loadtxt(path_16[i], delimiter='\t',
                            skiprows=2, usecols=(2))
    # q_16_x
    q_16 = psi_pa * np.loadtxt(path_16[i], delimiter='\t',
                               skiprows=2, usecols=(1))
    # Set color
    if (p_16[0] > (7000*psi_pa)):
        c = colors[4]
        marker = markers[4]
        name = 'pconf = 8000Psi'
    elif (p_16[0] > (4000*psi_pa)):
        c = colors[3]
        marker = markers[3]
        name = 'pconf = 5000Psi'
    elif (p_16[0] > (2300*psi_pa)):
        c = colors[2]
        marker = markers[2]
        name = 'pconf = 2500Psi'
    elif (p_16[0] > (800*psi_pa)):
        c = colors[1]
        marker = markers[1]
        name = 'pconf = 1000Psi'
    else:
        c = colors[0]
        marker = markers[0]
        name = 'pconf = 500Psi'
    # plot data
    plt.scatter(astrain_16, q_16, color=c, s=1, marker=marker, edgecolor=c,
                alpha=1, label='Traxial test data 16-'+str(i+1)+", "+name)
# Plot q-straina
for i in range(simulation_data):
    astrain_simulation = np.loadtxt(
        path_simulation_strain[i], delimiter='\t', usecols=(2))
    p_simulation = np.loadtxt(
        path_simulation_p_q[i], delimiter='\t', usecols=(0))
    q_simulation = np.loadtxt(
        path_simulation_p_q[i], delimiter='\t', usecols=(1))
    plt.plot(-astrain_simulation, q_simulation, 'k-',
             alpha=0.1)
    # Set color
    if (p_simulation[0] > (7000*psi_pa)):
        c = colors[4]
        name = 'pconf = 8000Psi'
    elif (p_simulation[0] > (4000*psi_pa)):
        c = colors[3]
        name = 'pconf = 5000Psi'
    elif (p_simulation[0] > (2300*psi_pa)):
        c = colors[2]
        name = 'pconf = 2500Psi'
    elif (p_simulation[0] > (800*psi_pa)):
        c = colors[1]
        name = 'pconf = 1000Psi'
    else:
        c = colors[0]
        name = 'pconf = 500Psi'
    plt.plot(-astrain_simulation, q_simulation, color=c, alpha=1.0, linewidth=3,
             label='BMCC Model'+", "+name)


# Legend
legend = plt.legend(prop=font_legend)
# Save figure
plt.savefig('q-straina.png')
######################################################
# e-astrain curves
# plot
plt.figure()
figsize = 11, 9
figure, ax = plt.subplots(figsize=figsize)
# Axial ranges
plt.axis([0, 0.06, 0.3, 0.7])
plt.tick_params(labelsize=16)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
# axis label name
plt.xlabel('axial strain', font_label)
plt.ylabel('void ratio', font_label)
# Compute e0
e0 = [0] * 5
e0[0] = data['e_ref'] - data['lambda'] * \
    math.log(8000*psi_pa / data['p_ref']) - \
    data['kappa'] * math.log(data['pc0']/(8000*psi_pa))
e0[1] = data['e_ref'] - data['lambda'] * \
    math.log(5000*psi_pa / data['p_ref']) - \
    data['kappa'] * math.log(data['pc0']/(5000*psi_pa))
e0[2] = data['e_ref'] - data['lambda'] * \
    math.log(2500*psi_pa / data['p_ref']) - \
    data['kappa'] * math.log(data['pc0']/(2500*psi_pa))
e0[3] = data['e_ref'] - data['lambda'] * \
    math.log(1000*psi_pa / data['p_ref']) - \
    data['kappa'] * math.log(data['pc0']/(1000*psi_pa))
e0[4] = data['e_ref'] - data['lambda'] * \
    math.log(500*psi_pa / data['p_ref']) - \
    data['kappa'] * math.log(data['pc0']/(500*psi_pa))
# Plot traxial test data
for i in range(triaxial_data):
    # p_16_x
    p_16 = psi_pa * np.loadtxt(path_16[i], delimiter='\t',
                               skiprows=2, usecols=(0))
    # astrain_16_x
    astrain_16 = np.loadtxt(path_16[i], delimiter='\t',
                            skiprows=2, usecols=(2))
    # rstrain1_16_x
    rstrain1_16 = np.loadtxt(path_16[i], delimiter='\t',
                             skiprows=2, usecols=(3))
    # rstrain2_16_x
    rstrain2_16 = np.loadtxt(path_16[i], delimiter='\t',
                             skiprows=2, usecols=(4))
    # Initialise e
    e_16 = [0] * len(astrain_16)
    # Set color
    if (p_16[0] > (7000*psi_pa)):
        c = colors[4]
        marker = markers[4]
        name = 'pconf = 8000Psi'
        e_16[0] = e0[0]
    elif (p_16[0] > (4000*psi_pa)):
        c = colors[3]
        marker = markers[3]
        name = 'pconf = 5000Psi'
        e_16[0] = e0[1]
    elif (p_16[0] > (2300*psi_pa)):
        c = colors[2]
        marker = markers[2]
        name = 'pconf = 2500Psi'
        e_16[0] = e0[2]
    elif (p_16[0] > (800*psi_pa)):
        c = colors[1]
        marker = markers[1]
        name = 'pconf = 1000Psi'
        e_16[0] = e0[3]
    else:
        c = colors[0]
        marker = markers[0]
        name = 'pconf = 500Psi'
        e_16[0] = e0[4]
    for j in range(len(astrain_16)-1):
        e_16[j+1] = e_16[j] - ((rstrain1_16[j+1] - rstrain1_16[j] + rstrain2_16[j+1] -
                                rstrain2_16[j] + astrain_16[j+1]-astrain_16[j]) * (1 + e_16[0]))

    # plot data
    plt.scatter(astrain_16, e_16, color=c, s=5, marker=marker,
                alpha=1, label='Traxial test data 16-'+str(i+1)+", "+name)

# Plot e-astrain
for i in range(simulation_data):
    p_simulation = np.loadtxt(
        path_simulation_p_q[i], delimiter='\t', usecols=(0))
    astrain_simulation = np.loadtxt(
        path_simulation_strain[i], delimiter='\t', usecols=(2))
    e_simulation = np.loadtxt(
        path_simulation_e[i], delimiter='\t', usecols=(0))
    # Set color
    if (p_simulation[0] > (7000*psi_pa)):
        c = colors[4]
        name = 'pconf = 8000Psi'
    elif (p_simulation[0] > (4000*psi_pa)):
        c = colors[3]
        name = 'pconf = 5000Psi'
    elif (p_simulation[0] > (2300*psi_pa)):
        c = colors[2]
        name = 'pconf = 2500Psi'
    elif (p_simulation[0] > (800*psi_pa)):
        c = colors[1]
        name = 'pconf = 1000Psi'
    else:
        c = colors[0]
        name = 'pconf = 500Psi'
    plt.plot(-astrain_simulation, e_simulation, color=c, alpha=1.0,
             label='BMCC Model'+", "+name)
# Legend
legend = plt.legend(prop=font_legend)
# Save figure
plt.savefig('e-straina.png')
######################################################
# e-p curves
# plot
plt.figure()
figsize = 11, 9
figure, ax = plt.subplots(figsize=figsize)
# Axial ranges
plt.axis([0, 100000000, 0.3, 0.7])
plt.tick_params(labelsize=16)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Times New Roman') for label in labels]
# axis label name
plt.xlabel('p\' / Pa', font_label)
plt.ylabel('void ratio', font_label)
# Plot traxial test data
for i in range(triaxial_data):
    # p_16_x
    p_16 = psi_pa * np.loadtxt(path_16[i], delimiter='\t',
                               skiprows=2, usecols=(0))
    # astrain_16_x
    astrain_16 = np.loadtxt(path_16[i], delimiter='\t',
                            skiprows=2, usecols=(2))
    # rstrain1_16_x
    rstrain1_16 = np.loadtxt(path_16[i], delimiter='\t',
                             skiprows=2, usecols=(3))
    # rstrain2_16_x
    rstrain2_16 = np.loadtxt(path_16[i], delimiter='\t',
                             skiprows=2, usecols=(4))
    # Initialise e
    e_16 = [0] * len(astrain_16)

    # Set color
    if (p_16[0] > (7000*psi_pa)):
        c = colors[4]
        marker = markers[4]
        name = 'pconf = 8000Psi'
        e_16[0] = e0[0]
    elif (p_16[0] > (4000*psi_pa)):
        c = colors[3]
        marker = markers[3]
        name = 'pconf = 5000Psi'
        e_16[0] = e0[1]
    elif (p_16[0] > (2300*psi_pa)):
        c = colors[2]
        marker = markers[2]
        name = 'pconf = 2500Psi'
        e_16[0] = e0[2]
    elif (p_16[0] > (800*psi_pa)):
        c = colors[1]
        marker = markers[1]
        name = 'pconf = 1000Psi'
        e_16[0] = e0[3]
    else:
        c = colors[0]
        marker = markers[0]
        name = 'pconf = 500Psi'
        e_16[0] = e0[4]

    for j in range(len(astrain_16)-1):
        e_16[j+1] = e_16[j] - ((rstrain1_16[j+1] - rstrain1_16[j] + rstrain2_16[j+1] -
                                rstrain2_16[j] + astrain_16[j+1]-astrain_16[j]) * (1 + e_16[0]))

    # plot data
    plt.scatter(p_16, e_16, color=c, s=5, marker=marker,
                alpha=1.0, label='Traxial test data 16-'+str(i+1)+", "+name)
# Plot e-astrain
for i in range(simulation_data):
    p_simulation = np.loadtxt(
        path_simulation_p_q[i], delimiter='\t', usecols=(0))
    e_simulation = np.loadtxt(
        path_simulation_e[i], delimiter='\t', usecols=(0))
    # Set color
    if (p_simulation[0] > (7000*psi_pa)):
        c = colors[4]
        name = 'pconf = 8000Psi'
    elif (p_simulation[0] > (4000*psi_pa)):
        c = colors[3]
        name = 'pconf = 5000Psi'
    elif (p_simulation[0] > (2300*psi_pa)):
        c = colors[2]
        name = 'pconf = 2500Psi'
    elif (p_simulation[0] > (800*psi_pa)):
        c = colors[1]
        name = 'pconf = 1000Psi'
    else:
        c = colors[0]
        name = 'pconf = 500Psi'
    plt.plot(p_simulation, e_simulation, ':', color=c, alpha=1.0, linewidth=2,
             label='BMCC Model'+", "+name)
# Compute critical line
e_critical = [0]*len(p_yield)


# pc0
pc0 = data['pc0']
# lambda
lam = data['lambda']
# p_ref
pref = data['p_ref']
# kappa
kappa = data['kappa']

for i in range(len(p_yield)):
    if (p_yield[i]) == 0:
        p_yield[i] = 1E-5
    # if(p_yield[i] > pc0):
    #    e_critical[i] = e0-lam * math.log(p_yield[i] / pc0)
    # else:
    #    e_critical[i] = e0-lam * math.log(p_yield[i] / pc0)
    e_critical[i] = 1-lam * math.log(p_yield[i] / pref)

plt.plot(p_yield, e_critical, color='k', alpha=1.0, linewidth=2,
         label='Critical line')

# Legend
legend = plt.legend(prop=font_legend)
# Save figure
plt.savefig('e-p.png')
