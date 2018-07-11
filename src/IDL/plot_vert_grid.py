import numpy as np
import matplotlib.pyplot as plt

n_level = 26
levels = np.array(range(n_level))
nodes = np.exp(levels * 0.2 - 5)
faces = np.zeros(n_level)
faces[0:n_level-1] = (nodes[0:n_level-1] + nodes[1:n_level]) / 2.
faces[n_level-1] = nodes[n_level-1] + (nodes[n_level-1] - nodes[n_level-2]) / 2.

sim_profiles = np.genfromtxt(
    './test_grid/sim_profiles.csv', dtype=None, skip_header=2, delimiter=',',
    names=['cos', 'co2'])

fig = plt.figure(figsize=(6, 5))
plt.subplot(121)
plt.plot([0, 0], [0, faces[-1]], linewidth=1, color='k')
plt.ylabel('Depth (m)')
plt.plot([2, 2], [0, faces[-1]], linewidth=1, color='k')
plt.axis([-0.5, 2.5, 1.2, -0.05])
plt.tick_params(axis='x', which='both', bottom='off', top='off',
                labelbottom='off')
plt.plot([0, 2], [0, 0], linewidth=1, color='k')
for i in levels:
    plt.plot([0, 2], [faces[i], faces[i]], linewidth=1, color='k')

plt.plot(np.ones(n_level), nodes, '*', color='b')

plt.subplot(122)
plt.plot([0, 0], [0, 3], linewidth=1, color='k')
plt.plot([2.5, 2.5], [0, 3], linewidth=1, color='k')
for i in range(4):
    plt.plot([0, 2.5], [i, i], linewidth=1, color='k')
plt.plot([1.25, 1.25, 1.25], [0.5, 1.5, 2.5], '*', color='b')
plt.axis([-0.25, 3.5, -1, 4])
plt.text(2.75, 0.5, '$i+1$', ha='left', va='center', fontsize=16)
plt.text(2.75, 1.5, '$i$', ha='left', va='center', fontsize=16)
plt.text(2.75, 2.5, '$i-1$', ha='left', va='center', fontsize=16)

plt.text(1.4, 0.4, '$C_{i+1}$', ha='left', va='center', fontsize=16)
plt.text(1.4, 1.4, '$C_i$', ha='left', va='center', fontsize=16)
plt.text(1.4, 2.4, '$C_{i-1}$', ha='left', va='center', fontsize=16)

plt.arrow(0.75, 1.2, 0, -0.4, head_width=0.1,
          head_length=0.1, fc='red', ec='red')
plt.arrow(0.75, 2.2, 0, -0.4, head_width=0.1,
          head_length=0.1, fc='red', ec='red')

plt.text(0.9, 1.8, r'$J_{i-1\rightarrow i}$',
         ha='left', va='center', color='red', fontsize=16)
plt.text(0.9, 0.8, r'$J_{i\rightarrow i+1}$',
         ha='left', va='center', color='red', fontsize=16)

plt.axis('off')

plt.text(-0.5, 3.75, '(a)', ha='center', va='center', fontsize=14)
plt.text(3.25, 3.75, '(b)', ha='center', va='center', fontsize=14)

# z_grid = np.insert(nodes, 0, 0.)
# plt.subplot(132)
# plt.plot(sim_profiles['cos'], z_grid, 'm*-')
# plt.ylabel('Depth (m)')
# plt.xlabel('COS (pptv)')
# plt.axis([0,550,1.2,-0.05])
# plt.plot([0,550], [0,0], 'k-', linewidth=1)
# plt.plot([0,550], [faces[12],faces[12]], 'k--')

# plt.subplot(133)
# plt.plot(sim_profiles['co2'], z_grid, 'r*-')
# plt.ylabel('Depth (m)')
# plt.xlabel('CO2 (ppmv)')
# plt.axis([0,5000,1.2,-0.05])
# plt.plot([0,5000], [0,0], 'k-', linewidth=1)
# plt.plot([0,5000], [faces[12],faces[12]], 'k--')

plt.tight_layout()
plt.savefig('./plots/mdl_plts/vert_grid.eps', dpi=300)
plt.savefig('./plots/mdl_plts/png/vert_grid.png', dpi=300)
# plt.show()
