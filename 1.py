# Import math Library
import math

from sympy import Identity, Matrix, MatrixSymbol, Symbol, init_printing, symbols

init_printing(use_unicode=True)

m= Matrix([1, 2, 3])
#print(m)

q1 = Symbol('q1')
q2 = Symbol('q2')
q3 = Symbol('q3')

X = MatrixSymbol('X', 3, 3)
Y = MatrixSymbol('Y', 3, 3)
(X.T*X).I*Y
print(X**(-1)*X.T**(-1)*Y)
print(Matrix(X))


# Return the sine value of 30 degrees
#print(round(math.cos(math.radians(-90))))


