$Title Stochastic fixed-charge transportation problem (SFCTP)

* relative optimality tolerance in solving MIP problems
option OptcR = 0, Decimals = 6

sets
   I         origins              / i1 * i4 /
   J         destinations         / j1 * j3 /
   S         scenarios            / s000 * s099 /

parameters
   A(i)      product offer        / i1 20, i2 30, i3 40, i4 20 /
   B(j)      product demand       / j1 21, j2 51, j3 31 /
   P(s)      scenario probability
   BS(s,  j) product demand stochastic ;

BS(s,j) = B(j) * [1+round(uniform(-0.05,0.05),2)] ;
P (s)   = 1/card(s) ;

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

loop (s, abort $(sum[i, A(i)] < sum[j, BS(s,j)]) 'Infeasible problem' )

positive variable
   X(s,i,j)     arc flow
binary   variable
   Y(i,j  )     arc investment decision
variables
   Z1           objective function

equations
   EQ_OBJ           complete problem objective function
   Offer    (s,i  ) offer  at origin
   Demand   (s,  j) demand at destination
   FlowLimit(s,i,j) arc flow limit ;

EQ_OBJ           .. Z1 =e= sum[(i,j), F(i,j)*Y(i,j)] + sum[(s,i,j), P(s)*C(i,j)*X(s,i,j)] ;
Offer    (s,i  ) .. sum[j, X(s,i,j)] =l= A (  i) ;
Demand   (s,  j) .. sum[i, X(s,i,j)] =g= BS(s,j) ;
FlowLimit(s,i,j) .. X(s,i,j) =l= 100 * Y(i,j) ;

model Complete      / EQ_OBJ Offer Demand FlowLimit / ;

X.up(s,i,j) = 100 ;

Complete.OptFile = 1 ;

file COPT / cplex.opt / ;
put  COPT putclose 'writelp FCTP_Sto.lp' / ;

solve Complete using MIP minimizing Z1

display Z1.l, Y.l
