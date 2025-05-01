import numpy as np

debug = False


def splitting(x0, xk, h):
    xs = []
    x = x0
    while x < xk:
        xs.append(x)
        x += h
    xs.append(xk)
    return xs


# Рунге-Кутт 4го порядка
p = 4
as_ = [0, 0.5, 0.5, 1]
bs = [[0.5], [0, 0.5], [0, 0, 1]]
cs = [1/6, 1/3, 1/3, 1/6]

def getKs(x: float, y: np.ndarray, h):
    dim = y.shape[0]
    Ks = np.empty((p, dim))

    for i in range(p):
        newX = x + as_[i] * h
        newY = np.copy(y)
        for j in range(i):
            newY += bs[i - 1][j] * Ks[j]

        K = h * f(newX, newY)
        if debug: print(f"\tK{i + 1} = {K}")
        Ks[i] = K

    return Ks

def getDeltaY(x: float, y: np.ndarray, h):
    Ks = getKs(x, y, h)
    dim = Ks.shape[1]
    sum_ = np.zeros(dim)
    for i in range(p):
        sum_ += cs[i] * Ks[i]
    if debug: print(f"\tdeltaY = {sum_}")
    return sum_

def RungeKutta(xs: list, y0: np.ndarray, h):
    N = len(xs) - 1
    dim = y0.shape[0]
    ys = np.empty((N + 1, dim))
    ys[0] = y0

    if debug: print(f"N = {N}, dim = {dim}")

    for k in range(1, N + 1):
        if debug: print(f"Шаг {k}")
        ys[k] = ys[k - 1] + getDeltaY(xs[k - 1], ys[k - 1], h)
        if debug: print(f"\ty = {ys[k]}")

    return ys


def Euler(xs: list, y0: np.ndarray, h):
    N = len(xs) - 1
    dim = y0.shape[0]
    ys = np.empty((N + 1, dim))
    ys[0] = y0

    for k in range(N):
        ys[k + 1] = ys[k] + h * f(xs[k], ys[k])

    return ys


def Adams(xs: list, y0s: np.ndarray, h):
    N = len(xs) - 1
    dim = y0s.shape[1]
    ys = np.empty((N + 1, dim))
    
    fs = np.empty((N + 1, dim))
    for i in range(4):
        ys[i] = np.copy(y0s[i])
        fs[i] = f(xs[i], ys[i])

    for k in range(4, N + 1):
        ys[k] = ys[k - 1] + h/24 * (55 * fs[k - 1] - 59 * fs[k - 2] + 37 * fs[k - 3] - 9 * fs[k - 4])
        fs[k] = f(xs[k], ys[k])

    return ys


def RungeError(ys: np.ndarray, ys2: np.ndarray, p):
    k = 2
    error = 0
    for i in range(ys.shape[0]):
        error = max(error, abs(ys2[i * 2][0] - ys[i][0]) / (k ** p - 1))
    return error


# Функции
def f(x: float, y: np.ndarray):
    return np.array([
        y[1],
        ((x + 1) * y[1] - y[0]) / x
    ])

def getTrueY(x):
    return x + 1 + np.exp(x)

# Интервал
a = 1
b = 2
# Начальные условия
y0 = np.array([2 + np.e, 1 + np.e])


# Считаем для шага из варианта

h = 0.1

xs = splitting(a, b, h)
ysEuler = Euler(xs, y0, h)
ysRungeKutta = RungeKutta(xs, y0, h)
ysAdams = Adams(xs, ysRungeKutta, h)

print(f"Шаг: {h}")
for i in range(len(xs)):
    y = getTrueY(xs[i])

    print(f"xk = {np.round(xs[i], 5)}, y(xk) = {np.round(y, 5)}")

    errorEuler = abs(ysEuler[i][0] - y)
    errorRungeKutta = abs(ysRungeKutta[i][0] - y)
    errorAdams = abs(ysAdams[i][0] - y)

    print(f"\tЭйлер:      yk = {np.round(ysEuler[i][0], 5)}, e = {np.round(errorEuler, 8)}")
    print(f"\tРунге-Кутт: yk = {np.round(ysRungeKutta[i][0], 5)}, e = {np.round(errorRungeKutta, 8)}")
    print(f"\tАдамс:      yk = {np.round(ysAdams[i][0], 5)}, e = {np.round(errorAdams, 8)}")


# Считаем для шага в два раза короче, чтобы применить оценку Рунге

h2 = h / 2

xs2 = splitting(a, b, h2)
ysEuler2 = Euler(xs2, y0, h2)
ysRungeKutta2 = RungeKutta(xs2, y0, h2)
ysAdams2 = Adams(xs2, ysRungeKutta2, h2)

print("===================================================================")
print("Апостериорные оценки погрешности по Рунге:")
print(f"\tЭйлер:      {RungeError(ysEuler, ysEuler2, 1)}")
print(f"\tРунге-Кутт: {RungeError(ysRungeKutta, ysRungeKutta2, 4)}")
print(f"\tАдамс:      {RungeError(ysAdams, ysAdams2, 3)}")
