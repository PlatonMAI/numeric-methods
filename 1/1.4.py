import numpy as np
import sys

def isSym(A):
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                print("Матрица не симметрична!")
                return False
    return True


def rotationMethod(A, eps):
    n = len(A)
    A = np.copy(A)
    U = np.eye(n)
    iter = 0
    while (True):
        iter += 1

        # Выбрали максимальный элемент
        I, J = -1, -1
        max_ = 0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > max_:
                    I, J = i, j
                    max_ = abs(A[i][j])
        
        # Строим Uk
        Uk = np.eye(n)
        phi = 0.5 * np.atan((2 * A[I][J]) / (A[I][I] - A[J][J]))
        Uk[I][I] = np.cos(phi)
        Uk[I][J] = -np.sin(phi)
        Uk[J][I] = np.sin(phi)
        Uk[J][J] = np.cos(phi)

        # Строим U
        U = np.dot(U, Uk)

        # Строим Ak
        A = np.dot(np.dot(np.transpose(Uk), A), Uk)

        # Проверяем критерий окончания
        t = 0
        for i in range(n):
            for j in range(i + 1, n):
                t += A[i][j] ** 2
        t = np.sqrt(t)
        if (t < eps):
            break

    lambdas = np.array([A[i][i] for i in range(n)])
    return lambdas, U, iter


n = int(input())

A = np.empty((n, n))
for i in range(n):
    A[i] = np.array(list(map(int, input().split(" "))))

if not isSym(A):
    sys.exit()

eps = float(input())


lambdas, U, iter = rotationMethod(A, eps)

print("Количество итераций: ", iter)
print()

w, v = np.linalg.eig(A)
print("Вычисленные собственные значения: ", lambdas)
print("Проверка собственных значений: ", w)
print()

print("Вычисленные собственные векторы:\n", U)
print("Проверка собственных значений:\n", v)
