# f(x) = 4^x - 5x - 2 = 0

# x = log4(5x + 2)

# 1.5
# [1, 2]

import math

def f(x):
    return 4 ** x - 5 * x - 2
def fDer(x):
    return math.log1p(3) * 4 ** x - 5
def fDer2(x):
    return math.log1p(3) ** 2 * 4 ** x
def newton(x0, eps):
    if (f(x0) * fDer2(x0) <= 0):
        raise Exception()

    xPrev = x0
    iter = 0
    while (True):
        iter += 1
        xCur = xPrev - f(xPrev) / fDer(xPrev)
        if abs(xCur - xPrev) < eps:
            break
        xPrev = xCur
    return xCur, iter


def phi(x):
    return math.log(5 * x + 2, 4)
def simpleIterations(x0, q, eps):
    xPrev = x0
    iter = 0
    while (True):
        iter += 1
        xCur = phi(xPrev)
        error = q / (1 - q) * abs(xCur - xPrev)
        if error < eps:
            break
        xPrev = xCur
    return xCur, iter


eps = float(input("Точность: "))

x0 = 0.5

newtonAns, iter = newton(x0, eps)
print("Метод Ньютона")
print("\tКорень: ", newtonAns)
print("\tКоличество итераций: ", iter)

q = 5 / (math.log1p(3) * 7)
simpleIterationsAns, iter = simpleIterations(x0, q, eps)
print("Метода простых итераций")
print("\tКорень: ", simpleIterationsAns)
print("\tКоличество итераций: ", iter)
