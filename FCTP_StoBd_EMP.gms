$Title Fixed-charge transportation problem (FCTP) solved by Benders decomposition

* relative optimality tolerance in solving MIP problems
option OptcR = 0, Decimals = 6

sets
   L         iterations        / l001 * l200 /
   LL(l)     iterations subset
   I         origins           / i1 * i4  /
   J         destinations      / j1 * j3  /
   S         scenarios         / s000 * s099 /

parameters
   A(i)      product offer        / i1 20, i2 30, i3 40, i4 20 /
   B(j)      product demand       / j1 21, j2 51, j3 31        /
   P(s)      scenario probability
   YS(s,i,j) arc investment decision
   XS(s,i,j) arc flow
   NS(s,  j) demand not served
   PI(s,i,j) dual variables of second stage constraints
   PTY       penalty for demand not served / 1000 / ;

table C(i,j) per unit variable transportation cost
       j1  j2  j3
   i1   1   2   3
   i2   3   2   1
   i3   2   3   4
   i4   4   3   2 ;

table F(i,j) fixed transportation cost
       j1  j2  j3
   i1  10  20  30
   i2  20  30  40
   i3  30  40  50
   i4  40  50  60 ;

parameter
   BS(s,  j) product demand ;

BS(s,j) = B(j) * [1+round(uniform(-0.05,0.05),2)] ;
P (s)   = 1/card(s) ;

abort $(sum[i, A(i)] < sum[j, B(j)]) 'Infeasible problem'

parameters
   BdTol            relative Benders tolerance / 1e-6 /
   Z_Lower          lower bound                / -inf /
   Z_Upper          upper bound                /  inf /
   Y_L    (l,  i,j) first stage variables values               in iteration l
   PI_L   (l,s,i,j) dual variables of second stage constraints in iteration l
   Delta  (l)       cut type (feasibility 0 optimality 1)      in iteration l
   Z2_L   (l)       subproblem objective function value        in iteration l

positive variable
   X(i,j)       arc flow
   N(  j)       demand not served

binary   variable
   Y(i,j)       arc investment decision

variables
   Z1           first  stage objective function
   Z2           second stage objective function
   Theta        recourse function

equations
   EQ_Z1          first  stage     objective function
   EQ_Z2          second stage     objective function
   EQ_OBJ         complete problem objective function
   Offer    (i  ) offer  at origin
   Demand   (  j) demand at destination
   FlowLimit(i,j) arc flow limit
   Bd_Cuts    (l) Benders cuts ;

EQ_Z1          .. Z1 =e= sum[(i,j), F(i,j)*Y(i,j)] + Theta ;
EQ_Z2          .. Z2 =e=                             sum[(i,j), C(i,j)*X(i,j)] + sum[j, PTY*N(j)] ;
EQ_OBJ         .. Z1 =e= sum[(i,j), F(i,j)*Y(i,j)] + sum[(i,j), C(i,j)*X(i,j)] + sum[j, PTY*N(j)] ;
Offer    (i  ) .. sum[j, X(i,j)]        =l= A(i) ;
Demand   (  j) .. sum[i, X(i,j)] + N(j) =g= B(j) ;
FlowLimit(i,j) .. X(i,j) =l= 100 * Y(i,j) ;
Bd_Cuts(ll)    .. Delta(ll) * Theta =g= Z2_L(ll) + sum[(s,i,j), (Y(i,j)-Y_L(ll,i,j)) * 100 * PI_L(ll,s,i,j)] ;

model Master_Bd     / EQ_Z1  Bd_Cuts                / ;
model Subproblem_Bd / EQ_Z2  Offer Demand FlowLimit / ;
model Complete      / EQ_OBJ Offer Demand FlowLimit / ;

X.up(i,j) = 100 ;

* EMP annotations
file emp / '%emp.info%' / ; emp.pc=2 ; emp.pw=1020
* define probability and values of the stochastic parameter
put  emp '* problem %GAMS.i%' / 'jrandvar '
loop (j,
   put B.tn(j) ' '
)
loop (s,
   put P(s)
   loop (j,
      put BS(s,j) /
   )
)
* define stochastic parameter, variable and constraints of the second stage
putclose emp / 'stage 2 B X Offer Demand FlowLimit'

set dict / s         . scenario . ''
           B         . randvar  . BS
           X         . level    . XS
           N         . level    . NS
           Y         . level    . YS
           FlowLimit . marginal . PI / ;

loop (s, abort $(sum[i, A(i)] < sum[j, BS(s,j)]) 'Infeasible problem' ) ;

* to allow CPLEX correctly detect rays in an infeasible problem
* only simplex method can be used and no preprocessing neither scaling options
* optimality and feasibility tolerances are very small to avoid primal degeneration

Subproblem_Bd.OptFile = 1 ;

file COPT / cplex.opt / ;
put  COPT putclose 'writelp FCTP_StoBd.lp' / 'names 1' / ;

file DOPT / de.opt / ;
put  DOPT putclose 'subsolver CPLEX' / 'subsolveropt 1' / ;

* parameter initialization

LL   (l)       = no ;
Delta(l)       =  0 ;
Z2_L (l)       =  0 ;
PI_L (l,s,i,j) =  0 ;
Y_L  (l,  i,j) =  0 ;

* Benders algorithm iterations
Theta.fx    =  0 ;
loop (l $(abs(1-Z_Lower/Z_Upper) > BdTol),

*  solving master problem
   solve Master_Bd using MIP minimizing Z1 ;

*  storing the master solution
   Y_L(l,i,j) = Y.l(i,j) ;

*  fixing first-stage variables and solving subproblem
   Y.fx( i,j) = Y.l(i,j) ;

*  solving subproblem

   solve Subproblem_Bd minimizing Z2 using emp scenario dict

*  storing parameters to build a new Benders cut
   if (Subproblem_Bd.ModelStat = 4,
      Delta(l) = 0 ;
      Z2_L (l) = Subproblem_Bd.SumInfes ;
   else
*     updating lower and upper bound
      Z_Lower =              Z1.l ;
      Z_Upper = min(Z_Upper, Z1.l - Theta.l + Z2.l) ;

      Theta.lo = -inf ;
      Theta.up =  inf ;

      Delta(l) = 1 ;
      Z2_L (l) = Z2.l ;
   ) ;

   PI_L(l,s,i,j) = PI(s,i,j) ;

   Y.lo(    i,j) = 0 ;
   Y.up(    i,j) = 1 ;

*  increase the set of Benders cuts
   LL(l) = yes ;
) ;

display Z_Lower, Z_Upper, Y.l
