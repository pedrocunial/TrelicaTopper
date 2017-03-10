"""Exercise 07/02."""
from math import sqrt

import numpy as np

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

def calc_degree():
    d = []
    for i in range(len(INCI)):
        x = [INCI[i][0] * 2 - 1, INCI[i][0] * 2,
             INCI[i][1] * 2 - 1, INCI[i][1] * 2]
        d.append(x)
    return d

degrees = calc_degree()


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
    K_e = []
    for i in range(len(INCI)):
        cos = INCI[i][3]
        sen = INCI[i][4]
        k_tmp = [
            [cos**2, cos*sen, -cos**2, -cos * sen],
            [cos*sen, sen**2, -cos * sen, -sen**2],
            [-cos**2, -cos * sen, cos**2, cos * sen],
            [-cos * sen, -sen**2, cos * sen, sen**2],
        ]
        for j in range(len(k_tmp)):
            k_tmp[j] = [(METER[i][0] * PROP[i][0] / INCI[i][2]) * k_tmp[j][w] for w in range(len(k_tmp[j]))]

        K_e.append(k_tmp)
    return K_e

calc_trelica()
K_e = calc_Ke()

def calc_kg():
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
    real_deal = []
    for i in range(len(CORDS) * 2):
        real_deal.append(np.zeros(len(CORDS) * 2))

    for k in K_g:
        for i in range(len(k)):
            for j in range(len(k[i])):
                real_deal[i][j] += k[i][j]

    return real_deal

real_deal = calc_real_deal()
print(real_deal)

P_g = [0, 0, 0, 0, 0, 0, 0, -1000]
