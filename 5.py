import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


x, y = sp.symbols("x y")

a = 0
b = 1.8
c = 1.515
k = (0.4, 1.4)
q = (3.2, 12)
f = 8 * x * (2 - x**2)
h = 0.1

def kx(x):
    return k[0] if x < c and x > a else k[1]

def qx(x):
    return q[0] if x < c and x > a else q[1]

fx = sp.lambdify(x, f, "numpy")

def TDMAsolver(a, b, c, d):
    nf = len(d) 
    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def get_diagonales(H):
    def xi(i):
        return a + i * H

    n = int((b - a) / H)

    aa = np.zeros(n-1, dtype=float)
    bb = np.zeros(n, dtype=float)
    cc = np.zeros(n-1, dtype=float)
    dd = np.zeros(n, dtype=float)

    bb[0] = kx(a) / H + 1/2
    cc[0] = -kx(a) / H
    dd[0] = 0
    for i in range(1, n-1):
        aa[i-1] = (1 / H**2) + qx(xi(i)) - (kx(xi(i)) / H)
        bb[i] = (kx(xi(i)) / H) - (2 / H**2)
        cc[i] = 1 / H**2
        dd[i] = fx(xi(i))
    aa[-1] = -kx(b) / H
    bb[-1] = kx(b) / H + 1/2
    dd[-1] = 0

    return aa, bb, cc, dd


def solve(H):
    a, b, c, d = get_diagonales(H)
    return TDMAsolver(a, b, c, d)


step = 1
while True: 
    res1 = solve(h)
    res2 = solve(h/2)
    dif = []
    for i in range(len(res1)):
        dif.append(abs(res1[i] - res2[i*2]))
    print(max(dif))
    if max(dif) < 0.001:
        print(res2, "\nStep=",step, sep='')
        h = h / 2
        break
    
    h = h / 2
    valuesx = np.arange(a, b, h)
    plt.plot(valuesx, res2, label=f'res h={h}')
    step = step+1
    

valuesx = np.arange(a, b, h)
plt.plot(valuesx, res2, label=f'res h={h}')
plt.legend()
plt.show()




