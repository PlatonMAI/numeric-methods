import numpy as np

NORM_ORD = np.inf

def buildMatrixAlphaAndBeta(A, b):
    n = len(b)
    alpha = np.empty((n, n))
    beta = np.empty(n)
    for i in range(n):
        beta[i] = b[i] / A[i][i]
        for j in range(n):
            alpha[i][j] = 0 if i == j else -A[i][j] / A[i][i]

    return alpha, beta

def FixedPointIteration(alpha, beta, eps):
    prevX = np.copy(beta)
    iter = 0
    norm = np.linalg.norm(alpha, ord=NORM_ORD)
    while (True):
        iter += 1
        curX = np.dot(alpha, prevX) + beta
        curEps = norm * (1 - norm) * np.linalg.norm(curX - prevX)
        prevX = curX
        if (curEps < eps):
            break

    return prevX, iter

def Seidel(A, b, eps):
    alpha, beta = buildMatrixAlphaAndBeta(A, b)

    B = np.zeros((n, n))
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i > j:
                B[i][j] = alpha[i][j]
            else:
                C[i][j] = alpha[i][j]

    tmp = np.linalg.inv((np.eye(n) - B))
    newAplha = np.dot(tmp, C)
    newBeta = np.dot(tmp, beta)

    return FixedPointIteration(newAplha, newBeta, eps)


n = int(input())

A = np.empty((n, n))
for i in range(n):
    A[i] = np.array(list(map(int, input().split(" "))))

b = np.array(list(map(int, input().split(" "))), dtype=float)

eps = float(input())

x, iter = Seidel(A, b, eps)

print("Число итераций: ", iter)
print("Решение системы: ", x)
print("Проверка решения: ", np.linalg.solve(A, b))
