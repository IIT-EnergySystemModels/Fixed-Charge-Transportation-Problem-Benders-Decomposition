sets
   I origins      / VIGO, ALGECIRAS /
   J destinations / MADRID, BARCELONA, VALENCIA / ;

parameters
   pA(i) origin capacity
       / VIGO      350
         ALGECIRAS 700 /

   pB(j) destination demand
       / MADRID    400
         BARCELONA 450
         VALENCIA  150 / ;

table pC(i,j) per unit transportation cost
          MADRID BARCELONA VALENCIA
VIGO       0.06     0.12     0.09
ALGECIRAS  0.05     0.15     0.11 ;

variables
   vX(i,j) units transported
   vCost   transportation cost

positive variable vX ;

equations
   eCost        transportation cost
   eCapacity(i) maximum capacity of each origin
   eDemand  (j) demand supply at destination ;

eCost        .. sum[(i,j), pC(i,j) * vX(i,j)] =e= vCost ;
eCapacity(i) .. sum[   j ,           vX(i,j)] =l= pA(i) ;
eDemand  (j) .. sum[ i   ,           vX(i,j)] =g= pB(j) ;

model mTransport / all / ;
solve mTransport using LP minimizing vCost
