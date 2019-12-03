import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

x, y = sp.symbols("x y")

E = 0.1
p = -1.
q = 2. * x**2
f = x + 1
a, b = 1.3, 2.4
c = [1., 0., 1.]
d = [1., 1., 1.]
h = 0.1

px = sp.lambdify(x, p, "numpy")
qx = sp.lambdify(x, q, "numpy")
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
        return i*H + a


    k = int((b-a)/H) + 1
    print(k)
    print(len(np.arange(a, b, H)))
    aa = np.zeros(k - 1, dtype=float)
    bb = np.zeros(k, dtype=float)
    cc = np.zeros(k - 1, dtype=float)
    dd = np.zeros(k, dtype=float)

    bb[0] = 1
    dd[0] = 1
    for i in range(1, k-1):
        aa[i-1] = (1 / H**2) + qx(xi(i)) - (px(xi(i)) / H)
        bb[i] = (px(xi(i)) / H) - (2 / H**2)
        cc[i] = 1 / H**2
        dd[i] = fx(xi(i))
    aa[-1] = 1
    bb[-1] = 1 + 1 / H
    dd[-1] = 3.2
    
    return aa, bb, cc, dd


def solve(H):
    A, B, C, D = get_diagonales(H)
    return TDMAsolver(A, B, C, D)


step = 1
while True: 
    res1 = solve(h)
    res2 = solve(h/2)
    dif = []
    for i in range(len(res1)):
        dif.append(abs(res1[i] - res2[i*2]))
    print(max(dif))
    if max(dif) < E:
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
