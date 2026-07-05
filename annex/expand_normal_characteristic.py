from sympy import (
    MatrixSymbol,
    BlockMatrix,
    symbols,
    block_collapse,
    pprint,
    init_printing,
    Eq,
    Rational,
    expand
)

init_printing(use_unicode=True)

dim = symbols("dim", integer=True, positive=True)


alpha = MatrixSymbol("alpha", 1, 1)
mu_e = MatrixSymbol("mu_e", dim, 1)
mu = BlockMatrix([[alpha], [mu_e]])

sigma2 = MatrixSymbol("sigma2", 1, 1)
gamma = MatrixSymbol("gamma", dim, 1)
p_e = MatrixSymbol("P_e", dim, dim)
p = BlockMatrix([[sigma2, gamma.T], [gamma, p_e]])

n = MatrixSymbol("n", 1, 1)
eta_e = MatrixSymbol("eta_e", dim, 1)
eta = BlockMatrix([[n], [eta_e]])

print("\nTerm multiplied by i\n")
expr1 = mu.T * eta
pprint(Eq(expr1, block_collapse(expr1)))

print("\nMatrix expansion\n")
expr2 = Rational(1,2) * eta.T * p * eta 
pprint(Eq(expr2, block_collapse(expr2)))

print("\nGVM Product expansion\n")
a = MatrixSymbol("A", dim, dim)
xi = MatrixSymbol("xi", dim, 1)
beta = MatrixSymbol("beta", dim, 1)
m = MatrixSymbol("m", 1, 1)

t1 = a.T * xi + m * beta
expr3 = Rational(1, 2) * t1.T * t1

pprint(expand(expr3))
