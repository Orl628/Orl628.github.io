using LinearCovarianceModels
sigma = toeplitz(3)
W = ml_degree_witness(sigma)
verify(W)
S = [  4/5  -9/5 -1/25
      -9/5 79/16 25/24
     -1/25 25/24 17/16]
critical_points(W,S)
mle(W,S)
critical_points(W,S,only_positive_definite=false)

exit()
julia
using LinearAlgebra
import HomotopyContinuation
import HomotopyContinuation: dim, parameters, soluions
const HC = HomotopyContinuation

# function that turns a symmetric matrix into a vector of its upper
# entries
sym_to_vec(S) = (n = size(S, 1); [S[i, j] for i = 1:n for j = i:n])

# set up variables

HC.@var theta[1:3]

sigma = [theta[1] theta[2] theta[3]
         theta[2] theta[1] theta[2]
         theta[3] theta[2] theta[1]]

m = 3  #number of vars
n = 3  #matrix dimension
N = binomial(n+1, 2)

HC.@var k[1:N] s[1:N]

# set up matrices of variables S and K

S = Matrix{eltype(s[1])}(undef, n, n)
K = Matrix{eltype(k[1])}(undef, n, n)
l = 1
for i = 1:n, j = i:n
  S[i, j] = S[j, i] = s[l]
  K[i, j] = K[j, i] = k[l]
  l += 1
end
S, K

# compute score equtions

L = -tr(K * sigma) + tr(S * K * sigma * K)
nablaL = HC.differentiate(L, theta)
KSigmaI = vec(K * sigma - Matrix(I, n, n))

size(nablaL), size(KSigmaI)

# create system of equations to be solved

system = HC.System([nablaL; KSigmaI], variables = [theta; k],
                   parameters = s)
size(system)

# find one solution to the system by exchanging variables and params

theta_0 = randn(ComplexF64, length(theta))
sigma_0 = HC.evaluate(sigma, theta => theta_0)
K_0 = inv(sigma_0)
x_0 = [theta_0; sym_to_vec(K_0)]
exprs = system(x_0, s[1:N])[1:end-length(K_0)]
s_0 = HC.find_start_pair(HC.System(exprs, s[1:N]))[1]

# find the rest of the solutions wrt s_0 with monodromy

W = HC.monodromy_solve(system, x_0, s_0)

HC.solutions(W)

# alternative: suppose we already knew that MLD = 3. Then may pass this
# argument.

W2 = HC.monodromy_solve(system, x_0, s_0;
                        target_solutions_count = 3)
                        
# use witness to compute a specific case.                        

S = [  4/5  -9/5 -1/25
      -9/5 79/16 25/24
     -1/25 25/24 17/16]
                        
result = HC.solve(
        system,
        HC.solutions(W);
        start_parameters = s_0,
        target_parameters = sym_to_vec(S),
    )    
HC.solutions(result)
