import numpy as np
import matplotlib.pyplot as plt

countPoints = int(1e3)

def f(x):
    return np.exp(x) + x

def M(x):
    return np.exp(x)

def lagrange(xs, x):
    n = len(xs)
    res = np.zeros(countPoints)

    for i in range(n):
        resCur = np.array([f(xs[i])] * countPoints)
        for j in range(n):
            if (i == j):
                continue
            resCur *= (x - xs[j]) / (xs[i] - xs[j])
        res += resCur

    return res

def newton(xs, x):
    n = len(xs)
    res = f(xs[0])

    polynom = np.ones(countPoints)
    diffsPrev = [f(x) for x in xs]
    for i in range(2, n + 1): # сколько аргументов у разделенной разности
        diffsCur = []
        polynom *= x - xs[i - 2]
        for j in range(n - i + 1): # с какого икса начинаем
            diffsCur.append((diffsPrev[j] - diffsPrev[j + 1]) / (xs[j] - xs[j + i - 1]))
        res += polynom * diffsCur[0]
        diffsPrev = np.copy(diffsCur)

    return res

xs = np.array([float(i) for i in input().split(" ")])

x = np.linspace(xs[0], xs[-1], countPoints)

fs = f(x)
lagranges = lagrange(xs, x)
newtons = newton(xs, x)

plt.plot(x, fs, linestyle='-', color=(1, 0, 0), label=f"exp(x) + x")
plt.plot(x, lagranges, linestyle='-', color=(0, 1, 0), label=f"Лагранж")
plt.plot(x, newtons, linestyle='-', color=(0, 0, 1), label=f"Ньютон")

plt.legend()

plt.grid(True)
plt.show()
