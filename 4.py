
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

x, y, c, c1, c2 = sp.symbols("x y c c1 c2")

f = sp.sqrt(x) + 4
fx = sp.lambdify(x, f, "numpy")
a = 0.5
ua = 2
b = 1.5
ub = 5
h = (b - a) / 150
valuesx = np.arange(a, b+h, h)
def TDMA(a, b, c, d):
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

# Пусть стержень состоит из 2-x материалов с различными свойствами:
def get_solve_1(k1, k2):
    M = (b - a) / 2 + a
    k = lambda x: (k1 if x <= M else k2)
    xi = lambda i: (a + i * h)
   
    N = int((b - a) / h) + 1

    aa = np.zeros(N - 1, dtype=float)
    bb = np.zeros(N, dtype=float)
    cc = np.zeros(N - 1, dtype=float)
    dd = np.zeros(N, dtype=float)

    bb[0] = 1
    dd[0] = ua
    for i in range(1, N-1):
        aa[i-1] = -k(xi(i)) / h ** 2
        bb[i] = 2 * k(xi(i)) / h ** 2
        cc[i] = -k(xi(i)) / h ** 2
        dd[i] = fx(xi(i))

    bb[-1] = 1
    dd[-1] = b

    return TDMA(aa, bb, cc, dd)


res1 = get_solve_1(1, 2)
res2 = get_solve_1(3, 1)
plt.plot(valuesx, res1, label='k1 >> k2')
plt.plot(valuesx, res2, label='k1 << k2')
plt.legend()
plt.show()

def get_solve_2(k1, k2, k3):
    M1 = (b - a) / 3 + a
    M2 = 2 * (b - a) / 3 + a

    def k(x):
        if x >= a and x < M1:
            return k1
        if x >= M1 and x < M2:
            return k2
        else:
            return k3

    xi = lambda i: (a + i * h)
   
    N = int((b - a) / h) + 1

    aa = np.zeros(N - 1, dtype=float)
    bb = np.zeros(N, dtype=float)
    cc = np.zeros(N - 1, dtype=float)
    dd = np.zeros(N, dtype=float)

    bb[0] = 1
    dd[0] = ua
    for i in range(1, N-1):
        aa[i-1] = -k(xi(i)) / h ** 2
        bb[i] = 2 * k(xi(i)) / h ** 2
        cc[i] = -k(xi(i)) / h ** 2
        dd[i] = fx(xi(i))

    bb[-1] = 1
    dd[-1] = b

    return TDMA(aa, bb, cc, dd)

res1 = get_solve_2(1, 5, 10)
res2 = get_solve_2(10, 5, 1)
res3 = get_solve_2(10, 20, 10)
res4 = get_solve_2(80, 4, 80)
plt.plot(valuesx, res1, label='k1<k2<k3')
plt.plot(valuesx, res2, label='k1>k2>k3')
plt.plot(valuesx, res3, label='k1>k2>k3')
plt.plot(valuesx, res4, label='k1>k2>k3')
plt.legend()
plt.show()
