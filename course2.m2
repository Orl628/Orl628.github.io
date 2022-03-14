restart

R = QQ[u_0, u_1, u_2, u_12, x_1, x_2]

f = (u_1 + u_12)/x_1 + u_12/(x_1 + x_2 + 2) -
    (u_2 + u_12)/(x_1 + 1) -
    (u_0 + u_1 + u_2 + u_12)/(x_1 + x_2 + 1)
    
g = (u_2 + u_12)/x_2 + u_12/(x_1 + x_2 + 2) -
    (u_1 + u_12)/(x_2 + 1) -
    (u_0 + u_1 + u_2 + u_12)/(x_1 + x_2 + 1)
    
--denoms = x_1 * x_2 * (x_1 + x_2 + 2) * (x_1 + 1) *
        (x_2 + 1) * (x_1 + x_2 + 1)
        
--f = f * denoms
--g = g * denoms

denom1 = denominator f
denom2 = denominator g

denoms = denom1 * denom2

f = f * denoms
g = g * denoms

I = ideal(f,g)

I = sub(I, R)

J = ideal(denoms)

L = saturate(I, J)

Lsolve = sub(L, {u_0 => 0, u_1 => 0, u_2 => 0,
             u_12 => 11})
             
(entries gens Lsolve)_0_2

degree Lsolve

loadPackage "NumericalAlgebraicGeometry"

Lsolve = sub(Lsolve, CC[x_1, x_2])

solveSystem (entries gens Lsolve)_0

L

L = sub(L, frac(QQ[u_0, u_1, u_2, u_12])[x_1,x_2])

degree L

(entries gens L)_0_0
(entries gens L)_0_1
(entries gens L)_0_2
