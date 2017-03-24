#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 07/02."""
from math import sqrt
from inverter import NumericMethods

import numpy as np


class bcolors:
    """Colors to print."""

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


CORDS = [
    [0, 0],
    [0, 21],
    [21, 0],
    [21, 21]
]

INCI = [
    [1, 2],
    [1, 3],
    [3, 4],
    [2, 4],
    [2, 3],
    [1, 4]
]

PROP = [
    [1],
    [1],
    [1],
    [1],
    [sqrt(2)],
    [sqrt(2)]
]

METER = [
    [float("21e5")],
    [float("21e5")],
    [float("21e5")],
    [float("21e5")],
    [float("21e5")],
    [float("21e5")]
]

BC_NODES = [
    [1, 1],
    [1, 2],
    [2, 1],
    [2, 2]
]


def calc_degree():
    """Calculate degrees of freedom for each element."""
    d = []
    for i in range(len(INCI)):
        x = [INCI[i][0] * 2 - 1, INCI[i][0] * 2,
             INCI[i][1] * 2 - 1, INCI[i][1] * 2]
        d.append(x)
    return d


degrees = calc_degree()
print(degrees)


def calc_comp(x1, x2, y1, y2):
    """Calculate size element."""
    L = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return L


def calc_cos(x1, x2, y1, y2):
    """Calculate cos of element angle."""
    L = calc_comp(x1, x2, y1, y2)
    cos = (x2 - x1) / L
    return cos


def calc_sin(x1, x2, y1, y2):
    """Calculate sin of element angle."""
    L = calc_comp(x1, x2, y1, y2)
    sin = (y2 - y1) / L
    return sin


def calc_trelica():
    """Main function to calculate trelica."""
    for i in INCI:
        comp = calc_comp(CORDS[(i[0] - 1)][0], CORDS[(i[1] - 1)][0],
                         CORDS[(i[0] - 1)][1], CORDS[(i[1] - 1)][1])
        cos = calc_cos(CORDS[(i[0] - 1)][0], CORDS[(i[1] - 1)][0],
                       CORDS[(i[0] - 1)][1], CORDS[(i[1] - 1)][1])
        sin = calc_sin(CORDS[(i[0] - 1)][0], CORDS[(i[1] - 1)][0],
                       CORDS[(i[0] - 1)][1], CORDS[(i[1] - 1)][1])
        i.append(comp)
        i.append(cos)
        i.append(sin)


def calc_Ke():
    """Calculate Ke."""
    K_e = []
    for i in range(len(INCI)):
        cos = INCI[i][3]
        sen = INCI[i][4]
        k_tmp = [
            [cos**2, cos * sen, -cos**2, -cos * sen],
            [cos * sen, sen**2, -cos * sen, -sen**2],
            [-cos**2, -cos * sen, cos**2, cos * sen],
            [-cos * sen, -sen**2, cos * sen, sen**2],
        ]
        for j in range(len(k_tmp)):
            k_tmp[j] = [(METER[i][0] * PROP[i][0] / INCI[i][2]) *
                        k_tmp[j][w] for w in range(len(k_tmp[j]))]

        K_e.append(k_tmp)
    return K_e


calc_trelica()
K_e = calc_Ke()


def calc_kg():
    """Calculate Kg."""
    m = []
    for j in range(len(K_e)):
        arr = []
        for i in range(len(CORDS) * 2):
            arr.append(np.zeros(len(CORDS) * 2))
        m.append(arr)

    for i in range(len(K_e)):
        for j in range(len(degrees[i])):
            for d in range(len(degrees[i])):
                m[i][degrees[i][j] - 1][degrees[i][d] - 1] = K_e[i][j][d]
    return m


K_g = calc_kg()


def calc_real_deal():
    """Calculate real_deal."""
    real_deal = []
    for i in range(len(CORDS) * 2):
        real_deal.append(np.zeros(len(CORDS) * 2))

    for k in K_g:
        for i in range(len(k)):
            for j in range(len(k[i])):
                real_deal[i][j] += k[i][j]

    return np.matrix(real_deal)


real_deal = calc_real_deal()


P_g = [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    -1000
]

deleted_members = []


def cut(matrix):
    """Cut matrix."""
    global P_g
    global deleted_members
    P_g = np.matrix(P_g)
    _bc_nodes = np.matrix(BC_NODES)
    _bc_nodes = (_bc_nodes[np.argsort(_bc_nodes.A[:, 0])])[::-1]
    for line in _bc_nodes:
        deletable = 2 * line[:, 0] + (line[:, 1] - 2) - 1
        matrix = np.delete(matrix, deletable, 0)
        matrix = np.delete(matrix, deletable, 1)
        P_g = np.delete(P_g, deletable, 1)
        deleted_members.append(deletable)
    P_g = np.transpose(P_g)
    return matrix


matrix = cut(real_deal)


def calc_u():
    """Calculate U."""
    _m = np.asarray(matrix)
    u, _, _ = NumericMethods.gauss_seidel(100, 0.005, _m, P_g)
    u2 = np.zeros(len(real_deal))
    i = 0
    j = 0
    for i in range(len(real_deal)):
        if i not in deleted_members:
            u2[i] = u[j]
            j += 1
    return u2


U = calc_u()


def calc_strain(_id):
    """Deformation."""
    L = INCI[_id][2]
    cos, sin = INCI[_id][3], INCI[_id][4]
    ele_1, ele_2 = INCI[_id][0] - 1, INCI[_id][1] - 1
    cs_arr = np.array([
        -cos,
        -sin,
        cos,
        sin
    ])
    _u = np.array([
        U[ele_1 * 2],
        U[ele_1 * 2 + 1],
        U[ele_2 * 2],
        U[ele_2 * 2 + 1]
    ])
    return (1 / L) * np.dot(cs_arr, _u)


def calc_stress(_id):
    """Tensao."""
    return METER[_id][0] * calc_strain(_id)


deleted_members_reaction = []
Ur = U


def calc_reaction(matrix):
    """Calculate reaction."""
    global deleted_members_reaction
    global Ur
    _bc_nodes = np.matrix(BC_NODES)
    _bc_nodes = (_bc_nodes[np.argsort(_bc_nodes.A[:, 0])])[::-1]

    for line in _bc_nodes:
        deletable = 2 * line[:, 0] + (line[:, 1] - 2) - 1
        matrix = np.delete(matrix, deletable, 1)
        Ur = np.delete(Ur, deletable, 0)
        deleted_members_reaction.append(deletable)

    tamanho = len(matrix)
    for i in range(len(matrix) - 1, 0, -1):
        if i not in deleted_members_reaction:
            matrix = np.delete(matrix, i, 0)

    reac = np.dot(matrix, Ur)

    ur2 = np.zeros(tamanho)
    i = 0
    j = 0

    for i in range(tamanho):

        if i in deleted_members_reaction:
            ur2[i] = reac[::, j]
            j += 1

    return ur2


reacoes = calc_reaction(real_deal)


for i in range(len(METER)):
    print(str(bcolors.BOLD) + bcolors.OKBLUE +
          "=====================================")
    print(bcolors.HEADER + "Bar {}".format(i) + bcolors.ENDC)
    print(bcolors.WARNING +
          "{}Strain: {}{}".format(bcolors.BOLD, bcolors.ENDC, calc_strain(i)))
    print(bcolors.WARNING +
          "{}Stress: {}{}".format(bcolors.BOLD, bcolors.ENDC, calc_stress(i)))

print(str(bcolors.BOLD) + bcolors.OKBLUE +
      "=====================================")
print(bcolors.HEADER + "Reactions" + bcolors.ENDC)

# print(reacoes)

for i in range(len(reacoes)):
    if reacoes[i] != 0:
        print("{}{}R{}: {}{}".format(bcolors.WARNING,
                                     bcolors.BOLD, i,
                                     bcolors.ENDC, reacoes[i]))

n_i = 1
for i in range(len(U)):
    print(str(bcolors.BOLD) + bcolors.OKBLUE +
          "=====================================")
    print(bcolors.HEADER + "Node {}".format(i // 2) + bcolors.ENDC)
    print(bcolors.WARNING +
          "{}Displacement: {}{} in {}".format(bcolors.BOLD,
                                              bcolors.ENDC,
                                              U[i],
                                              "y" if i % 2 else "x"))
