# Deterministic Fixed-charge Transportation Problem (FCTP) solved by Benders decomposition

# Developed by

#    Andres Ramos
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingeniería - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu

#    MIT Energy Initiative
#    Massachusetts Institute of Technology
#    arght@mit.edu

#    May 2, 2017

# Define the packages
using JuMP    # used for mathematical programming
using Gurobi  # used as the solver

# system dimensions
L = 10  # iterations
I =  4  # origins
J =  3  # destinations
l =  1  # initial iteration

#  parameters
const A = [10, 30, 40, 20]                                  # product offer
const B = [20, 50, 30    ]                                  # product demand
const C = [ 1   2   3;  3   2   1;  2   3   4;  4   3   2]  # per unit variable transportation cost
const F = [10  20  30; 20  30  40; 30  40  50; 40  50  60]  #          fixed    transportation cost
const pDeficitCost = ones(J) * 1e4                          # penalty for the demand not supplied

if sum(A) < sum(B)
  @sprintf "Infeasible problem"
end

const BdTol =   1e-6  # relative Benders tolerance

Y_L   = zeros(L,I,J)  # first stage variables values               in iteration l
Y_J   = zeros(  I,J)  # first stage variables values               in iteration j
PI_L  = zeros(L,I,J)  # dual variables of second stage constraints in iteration l
Delta = zeros(L)      # cut type (feasibility 0 optimality 1)      in iteration l
Z2_L  = zeros(L)      # subproblem objective function value        in iteration l

Z_Lower = - 1e20      # lower bound (-Inf not allowed, fails in the if)
Z_Upper =   1e20      # upper bound ( Inf not allowed, fails in the if)

# fixed-charge transportation problem
Master_Bd     = Model(solver=GurobiSolver())
Subproblem_Bd = Model(solver=GurobiSolver())

# decision variables
@variable(Master_Bd,    - Inf <= Theta                 <= Inf           )  # recourse function
@variable(Master_Bd,             Y[       i=1:I,j=1:J], Bin             )  # arc investment decision
@variable(Subproblem_Bd,    0 <= X[       i=1:I,j=1:J] <= min(A[i],B[j]))  # arc flow
@variable(Subproblem_Bd,    0 <= vDeficit[      j=1:J] <=          B[j] )  # deficit of demand

# first and second stage objective function
@objective(Master_Bd,     :Min, sum(F[i,j] * Y[i,j] for i=1:I,j=1:J) + Theta                                                                              )
@objective(Subproblem_Bd, :Min,                                        sum(C[i,j] * X[i,j] for i=1:I,j=1:J) + sum(pDeficitCost[j] * vDeficit[j] for j=1:J))

@constraint(Subproblem_Bd, Offer[    i=1:I      ], sum(X[i,j] for j=1:J)               <=     A[i]               )  # offer  at origin
@constraint(Subproblem_Bd, Demand[         j=1:J], sum(X[i,j] for i=1:I) + vDeficit[j] >=          B[j]          )  # demand at destination
@constraint(Subproblem_Bd, FlowLimit[i=1:I,j=1:J],     X[i,j]                          <= min(A[i],B[j])*Y_J[i,j])  # arc flow limit

# fix initial value
setlowerbound(Theta, 0)
setupperbound(Theta, 0)

# Benders algorithm iterations
for l in 1:L
  if abs(1-Z_Lower/Z_Upper) > BdTol

    # solving master problem
    #print("\nThe current master problem model is \n", Master_Bd)
    Master_Bd_status = solve(Master_Bd)
    Z1   = getobjectivevalue(Master_Bd)

    for i in 1:I
      for j in 1:J
        # storing the master solution
        Y_L[l,i,j] = getvalue(Y[i,j])
        Y_J[  i,j] = getvalue(Y[i,j])
        chgConstrRHS(FlowLimit[i,j], min(A[i],B[j])*Y_J[i,j])
      end
    end

    # solving subproblem
    #print("\nThe current subproblem model is \n", Subproblem_Bd)
    Subproblem_Bd_status = solve(Subproblem_Bd)
    Z2       = getobjectivevalue(Subproblem_Bd)
    Z2_L[l]  = Z2

    # storing parameters to build a new Benders cut
    if Subproblem_Bd_status == :Infeasible
      # the problem has to be feasible because I am not able to obtain the sum of infeasibilities of the phase I
      Delta[l] =  0
    else
      # updating lower and upper bound
      Z_Lower =              Z1
      Z_Upper = min(Z_Upper, Z1 - getvalue(Theta) + Z2)
      println("Iteration ", l, " Status of the master problem is ", Master_Bd_status, " with o.f. Master_Bd ",     Z1, " Z_Lower = ", Z_Lower)
      println("Iteration ", l, " Status of the subproblem is ", Subproblem_Bd_status, " with o.f. Subproblem_Bd ", Z2, " Z_Upper = ", Z_Upper)

      setlowerbound(Theta, - Inf)
      setupperbound(Theta,   Inf)

      Delta[l] =  1
    end

    for i in 1:I
      for j in 1:J
        # storing the master solution
        PI_L[l,i,j] = getdual(FlowLimit[i,j])
      end
    end

    # Benders cuts
    @constraint(Master_Bd, Bd_Cuts[l], Delta[l] * Theta >= Z2_L[l] - sum(PI_L[l,i,j] * min(A[i],B[j]) * (Y_L[l,i,j] - Y[i,j]) for i=1:I,j=1:J))

  end
end

# complete problem
Complete      = Model(solver=GurobiSolver())

@variable(       Complete, 0 <= X[i=1:I,j=1:J] <= min(A[i],B[j]))
@variable(       Complete,      Y[i=1:I,j=1:J], Bin             )

@objective( Complete, :Min, sum(F[i,j] * Y[i,j] for i=1:I,j=1:J) + sum(C[i,j] * X[i,j] for i=1:I,j=1:J))

@constraint(Complete, Offer[   i=1:I       ], sum(X[i,j] for j=1:J) <=     A[i]             )
@constraint(Complete, Demand[         j=1:J], sum(X[i,j] for i=1:I) >=          B[j]        )
@constraint(Complete, FlowLimit[i=1:I,j=1:J],     X[i,j]            <= min(A[i],B[j])*Y[i,j])

#print("\nThe current complete problem model is \n", Complete)
Complete_status = solve(Complete)
Z1  = getobjectivevalue(Complete)
println("Status of the complete problem is ", Complete_status, " with o.f. Complete ", Z1)
