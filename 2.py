import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

x, y, c, c1, c2 = sp.symbols("x y c c1 c2")

p = sp.sin(2 * x)
q = 8 * (1 + sp.sin(x) ** 2)
f = 10 * sp.cos(x)
a, b = 1., 3.
UA, UB = 0, 0
E = 0.05
h = 0.1

px = sp.lambdify(x, p, "numpy")
qx = sp.lambdify(x, q, "numpy")
fx = sp.lambdify(x, f, "numpy")

def get_matrixA(H):
    def xi(i):
        return a + H * i

    k = int((b-a)/H) + 1

    res = np.zeros((k, k), dtype=float)

    res[0, 0] = 1
    for i in range(1, k-1):
        res[i, i-1] = (1 / H**2) + qx(xi(i)) - (px(xi(i)) / H)
        res[i, i] = (px(xi(i)) / H) - (2 / H**2)
        res[i, i+1] = 1 / H**2
    res[-1, -1] = 1

    return res

get_matrixA(0.5)
    

def get_matrixB(H):
    n = int((b-a) / H) + 1
    
    res = np.zeros(n, dtype=float)
    
    res[0] = a
    for i in range(1, int(n)-1):
        res[i] = fx(H*i+a)
    res[-1] = b
    
    return res

def solve(H):
    A = get_matrixA(H)
    B = get_matrixB(H)
    res = np.linalg.solve(A, B)
    
    return res

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
    step = step+1
    valuesx = np.arange(a, b+h, h)
    plt.plot(valuesx, res2, label=f'res h={h}')
    
valuesx = np.arange(a, b+h, h)
plt.plot(valuesx, res2, label=f'res h={h}')
plt.legend()
plt.show()
