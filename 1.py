import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

x, y, c, c1, c2 = sp.symbols("x y c c1 c2")
kx = x
fx = sp.sqrt(x)+4
a = 0.5
ua = 0
b = 1.5
ub = 5

C = (1, 2, 0.1, 1, 1, 1, 1)
Kx = (kx, c*kx, c*kx, 1/kx, kx, kx, kx)
UA = (ua, ua, ua, ua, -ua, ua, -ua)
UB = (ub, ub, ub, ub, ub, -ub, ub)

def Integrate(k, F):
    f1 = sp.integrate(F, x)
    f2 = -sp.integrate(f1 / k, x) + sp.integrate(c1 / k, x) + c2
    return f2
    
def solve1(F, a, ua, b, ub, C):
    f1 = sp.Eq(F.subs([(x, a), (c, C)]), ua)
    f2 = sp.Eq(F.subs([(x, b), (c, C)]), ub)
    return sp.solve([f1, f2])
    
res = Integrate(Kx[0], fx)
m = solve1(res, a, UA[0], b, UB[0], 1)
print("Решение уравнения ", res.subs([(c1, m[c1]), (c2, m[c2])]))

# Изменяя значения параметра с в коэффициенте теплопроводности,
# найти решения задачи для наборов параметров 1-3
solvs = []
for i in range(3):
    res = Integrate(Kx[i], fx)
    m = solve1(res, a, UA[i], b, UB[i], C[i])
    solvs.append(res.subs([(c1, m[c1]), (c2, m[c2]), (c, C[i])]))
    print("Для набора ", i + 1, " решение:\n", solvs[i])

# На одном чертеже построить графики найденных решений. Cравнить
# полученные результаты.

valuesx = np.arange(a, b, 0.01)

for i in range(3):
    np_solve = sp.lambdify(x, solvs[i], "numpy")
    plt.plot(valuesx, np_solve(valuesx), label = i)
plt.legend()
plt.show()

# Аналогично п.2, найти аналитическое решение для набора параметров 4.
# На одном чертеже построить графики решений для наборов 1 и 4. Cравнить
# полученные результаты.

res = Integrate(Kx[3], fx)
m = solve1(res, a, UA[3], b, UB[3], C[3])
solvs.append(res.subs([(c1, m[c1]), (c2, m[c2]), (c, C[3])]))
print("Решение уравнения 4: ", solvs[3])

np_solve1 = sp.lambdify(x, solvs[0], "numpy")
np_solve4 = sp.lambdify(x, solvs[3], "numpy")

plt.plot(valuesx, np_solve1(valuesx), label = 1)
plt.plot(valuesx, np_solve4(valuesx), label = 4)
plt.legend()
plt.show()

for i in range(3):
    res = Integrate(Kx[i+4], fx)
    m = solve1(res, a, UA[i+4], b, UB[i+4], C[i+4])
    solve = res.subs([(c1, m[c1]), (c2, m[c2]), (c, C[i])])
    print("Для набора ", i + 5, " решение:\n", solve)
    np_solve = sp.lambdify(x, solve, "numpy")
    plt.plot(valuesx, np_solve(valuesx), label = i + 5)
    
plt.legend()    
plt.show()
# Прим. графики 6 и 7 одинаковы, так как ua = 0