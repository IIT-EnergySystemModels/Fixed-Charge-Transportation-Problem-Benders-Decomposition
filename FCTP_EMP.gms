$title Deterministic fixed-charge transportation problem (FCTP)

$OnText
Developed by

   Andres Ramos
   Instituto de Investigacion Tecnologica
   Escuela Tecnica Superior de Ingenieria - ICAI
   UNIVERSIDAD PONTIFICIA COMILLAS
   Alberto Aguilera 23
   28015 Madrid, Spain
   Andres.Ramos@comillas.edu
   https://pascua.iit.comillas.edu/aramos/Ramos_CV.htm

   May 8, 2023

$OffText

* relative optimality tolerance in solving MIP problems
option OptcR = 0, Decimals = 6

sets
   I         origins           / i1 * i4 /
   J         destinations      / j1 * j3 / ;

parameters
   A(i)      product offer  / i1 20, i2 30, i3 40, i4 20 /
   B(j)      product demand / j1 21, j2 51, j3 31 / ;

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

positive variable
   X(i,j)       arc flow

binary   variable
   Y(i,j)       arc investment decision

variables
   Z1           objective function

equations
   EQ_OBJ         complete problem objective function
   Offer    (i  ) offer  at origin
   Demand   (  j) demand at destination
   FlowLimit(i,j) arc flow limit ;

EQ_OBJ         .. Z1 =e= sum[(i,j), F(i,j)*Y(i,j)] + sum[(i,j), C(i,j)*X(i,j)] ;
Offer    (i  ) .. sum[j, X(i,j)] =l= A(i) ;
Demand   (  j) .. sum[i, X(i,j)] =g= B(j) ;
FlowLimit(i,j) .. X(i,j) =l= 100 * Y(i,j) ;

model Complete / all / ;

X.up(i,j) = 100 ;

set S scenarios / s000 * s099 /
parameter
   BS(s,  j) product demand
   YS(s,i,j) arc investment decision
   XS(s,i,j) arc flow
   P (s    ) scenario probability ;

BS(s,j) = B(j) * [1+round(uniform(-0.05,0.05),2)] ;
P (s)   = 1/card(s) ;

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

set dict / s . scenario . ''
           B . randvar  . BS
           X . level    . XS
           Y . level    . YS / ;

loop (s, abort $(sum[i, A(i)] < sum[j, BS(s,j)]) 'Infeasible problem') ;

file COPT / cplex.opt / ;
put  COPT putclose 'writelp FCTP_EMP.lp' / 'names 1' / ;

file DOPT / de.opt / ;
put  DOPT putclose 'subsolver CPLEX' / 'subsolveropt 1' / ;

Complete.OptFile = 1 ;

solve Complete minimizing Z1 using emp scenario dict

display Z1.l, YS
