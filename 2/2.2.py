# 3x1 - cos(x2) = 0
# 3x2 - exp(x1) = 0

import numpy as np

NORM = np.inf

def getLU(A):
    n = len(A)
    LU = np.copy(A)
    swaps = []
    for k in range(n): # обнуляемый столбец
        if (LU[k][k] == 0): # Ищем ненулевой элемент
            ind = -1
            for i in range(k + 1, n):
                if LU[i][i] != 0:
                    ind = i
                    break
            if ind == -1:
                continue
            LU[[k, ind]] = LU[[ind, k]] # Меняем местами строки если нашли строку с ненулевым элементом
            swaps.append((k, ind))

        for i in range(k + 1, n): # текущая строка
            mu = LU[i][k] / LU[k][k]
            for j in range(k, n): # текущая столбец
                if (j == k):
                    LU[i][j] = mu
                else:
                    LU[i][j] -= mu * LU[k][j]

    return (LU, swaps)

def solverLU(LU, swaps, b):
    n = len(LU)

    b = np.copy(b)

    # Меняем строки в столбце свободных членов в соответствие с заменами строк в исходной матрице
    for swap in swaps:
        b[[swap[0], swap[1]]] = b[[swap[1], swap[0]]]

    # Lz = b
    z = np.zeros(n)
    for i in range(n):
        sum_ = sum([LU[i][j] * z[j] for j in range(i)])
        z[i] = b[i] - sum_

    # Ux = z
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        sum_ = sum([LU[i][j] * x[j] for j in range(n - 1, i, -1)])
        x[i] = (z[i] - sum_) / LU[i][i]

    return x


def f(x):
    return np.array([
        3 * x[0] - np.cos(x[1]),
        3 * x[1] - np.exp(x[0])
    ])
def fDer(x):
    return np.array([
        [3, np.sin(x[1])],
        [-np.exp(x[0]), 3]
    ])
def newton(x0, eps):
    xPrev = x0
    iter = 0
    while (True):
        iter += 1
        (LU, swaps) = getLU(fDer(xPrev))
        xDelta = solverLU(LU, swaps, -f(xPrev))
        xCur = xPrev + xDelta
        if np.linalg.norm(xCur - xPrev, NORM) < eps:
            break
        xPrev = xCur
    return xCur, iter


def phi(x):
    return np.array([
        np.cos(x[1]) / 3,
        np.exp(x[0]) / 3
    ])
def simpleIterations(x0, q, eps):
    xPrev = x0
    iter = 0
    while (True):
        iter += 1
        xCur = phi(xPrev)
        error = q / (1 - q) * np.linalg.norm(xCur - xPrev, NORM)
        if error < eps:
            break
        xPrev = xCur
    return xCur, iter


eps = float(input("Точность: "))

x0 = np.array([0, 0.3])

newtonAns, iter = newton(x0, eps)
print("Метод Ньютона")
print("\tКорень: ", newtonAns)
print("\tКоличество итераций: ", iter)

q = np.e / 3
simpleIterationsAns, iter = simpleIterations(x0, q, eps)
print("Метода простых итераций")
print("\tКорень: ", simpleIterationsAns)
print("\tКоличество итераций: ", iter)
