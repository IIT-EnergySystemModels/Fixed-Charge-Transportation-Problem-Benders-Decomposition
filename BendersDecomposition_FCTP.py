import pandas as pd
import pyomo.environ as pyo
from   pyomo.environ import ConcreteModel, Set, Param, Var, Binary, NonNegativeReals, RealSet, Constraint, Objective, minimize, Suffix, TerminationCondition
from   pyomo.opt     import SolverFactory

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

mFCTP      = ConcreteModel('Fixed-Charge Transportation Problem')
mMaster_Bd = ConcreteModel('Master problem')

mFCTP.i       = Set(initialize=['i1', 'i2', 'i3', 'i4'], doc='origins'     )
mFCTP.j       = Set(initialize=['j1', 'j2', 'j3'      ], doc='destinations')

mMaster_Bd.l  = Set(initialize=['it1', 'it2', 'it3', 'it4', 'it5', 'it6', 'it7', 'it8', 'it9', 'it10'], ordered=True, doc='iterations')
mMaster_Bd.ll = Set(                                                                                                  doc='iterations')

mFCTP.pA      = Param(mFCTP.i, initialize={'i1': 20, 'i2': 30, 'i3': 40, 'i4': 20}, doc='origin capacity'   )
mFCTP.pB      = Param(mFCTP.j, initialize={'j1': 20, 'j2': 50, 'j3':30           }, doc='destination demand')

FixedCost = {
    ('i1', 'j1'): 10,
    ('i1', 'j2'): 20,
    ('i1', 'j3'): 30,
    ('i2', 'j1'): 20,
    ('i2', 'j2'): 30,
    ('i2', 'j3'): 40,
    ('i3', 'j1'): 30,
    ('i3', 'j2'): 40,
    ('i3', 'j3'): 50,
    ('i4', 'j1'): 40,
    ('i4', 'j2'): 50,
    ('i4', 'j3'): 60,
    }

TransportationCost = {
    ('i1', 'j1'): 1,
    ('i1', 'j2'): 2,
    ('i1', 'j3'): 3,
    ('i2', 'j1'): 3,
    ('i2', 'j2'): 2,
    ('i2', 'j3'): 1,
    ('i3', 'j1'): 2,
    ('i3', 'j2'): 3,
    ('i3', 'j3'): 4,
    ('i4', 'j1'): 4,
    ('i4', 'j2'): 3,
    ('i4', 'j3'): 2,
    }

mFCTP.pF      = Param(mFCTP.i, mFCTP.j, initialize=FixedCost,          doc='fixed investment cost'       )
mFCTP.pC      = Param(mFCTP.i, mFCTP.j, initialize=TransportationCost, doc='per unit transportation cost')

mFCTP.vY      = Var  (mFCTP.i, mFCTP.j, bounds=(0,1), doc='units transported', within=Binary)
mMaster_Bd.vY = Var  (mFCTP.i, mFCTP.j, bounds=(0,1), doc='units transported', within=Binary)

mMaster_Bd.vTheta = Var(doc='transportation cost', within=RealSet)

mFCTP.vX      = Var  (mFCTP.i, mFCTP.j, bounds=(0.0,None), doc='units transported', within=NonNegativeReals)
mFCTP.vDNS    = Var  (         mFCTP.j, bounds=(0.0,None), doc='demand not served', within=NonNegativeReals)

def eCostMst(mMaster_Bd):
    return sum(mFCTP.pF[i,j]*mMaster_Bd.vY[i,j] for i,j in mFCTP.i*mFCTP.j) + mMaster_Bd.vTheta
mMaster_Bd.eCostMst = Objective(rule=eCostMst, sense=minimize, doc='total cost')

def eBd_Cuts(mMaster_Bd, ll):
    return mMaster_Bd.vTheta - Z2_L[ll] >= - sum(PI_L[ll,i,j] * min(mFCTP.pA[i],mFCTP.pB[j]) * (Y_L[ll,i,j] - mMaster_Bd.vY[i,j]) for i,j in mFCTP.i*mFCTP.j)

def eCostSubp(mFCTP):
    return sum(mFCTP.pC[i,j]*mFCTP.vX[i,j] for i,j in mFCTP.i*mFCTP.j) + sum(mFCTP.vDNS[j]*1000 for j in mFCTP.j)
mFCTP.eCostSubp = Objective(rule=eCostSubp, sense=minimize, doc='transportation cost')

def eCapacity(mFCTP, i):
    return sum(mFCTP.vX[i,j] for j in mFCTP.j) <= mFCTP.pA[i]
mFCTP.eCapacity  = Constraint(mFCTP.i,        rule=eCapacity,  doc='maximum capacity of each origin')

def eDemand  (mFCTP, j):
    return sum(mFCTP.vX[i,j] for i in mFCTP.i) + mFCTP.vDNS[j] >= mFCTP.pB[j]
mFCTP.eDemand    = Constraint(        mFCTP.j, rule=eDemand,    doc='demand supply at destination'   )

def eFlowLimit(mFCTP, i, j):
    return mFCTP.vX[i,j] <= min(mFCTP.pA[i],mFCTP.pB[j])*mFCTP.vY[i,j]
mFCTP.eFlowLimit = Constraint(mFCTP.i*mFCTP.j, rule=eFlowLimit, doc='arc flow limit'                 )

mFCTP.dual = Suffix(direction=Suffix.IMPORT)

Solver = SolverFactory('gurobi')
Solver.options['LogFile'] = 'mFCTP.log'

# initialization
Z_Lower = float('-inf')
Z_Upper = float(' inf')
BdTol   = 1e-6

Y_L   = pd.Series([0.]*len(mMaster_Bd.l*mFCTP.i*mFCTP.j), index=pd.MultiIndex.from_tuples(mMaster_Bd.l*mFCTP.i*mFCTP.j))
PI_L  = pd.Series([0.]*len(mMaster_Bd.l*mFCTP.i*mFCTP.j), index=pd.MultiIndex.from_tuples(mMaster_Bd.l*mFCTP.i*mFCTP.j))
Z2_L  = pd.Series([0.]*len(mMaster_Bd.l                ), index=mMaster_Bd.l)
Delta = pd.Series([0.]*len(mMaster_Bd.l                ), index=mMaster_Bd.l)

# Benders algorithm
mMaster_Bd.vTheta.fix(0)
for l in mMaster_Bd.l:
    if abs(1-Z_Lower/Z_Upper) > BdTol or l == mMaster_Bd.l.first():

        # solving master problem
        SolverResultsMst = Solver.solve(mMaster_Bd)
        Z1               = mMaster_Bd.eCostMst.expr()

        for i,j in mFCTP.i*mFCTP.j:
            # storing the master solution
            Y_L[l,i,j] = pyo.value(mMaster_Bd.vY[i,j])
            # fix investment decision for the subproblem
            mFCTP.vY[i,j].fix(Y_L[l,i,j])

        # solving subproblem
        SolverResultsSbp    = Solver.solve(mFCTP)
        Z2                  = mFCTP.eCostSubp.expr()
        Z2_L[l]             = Z2

        # storing parameters to build a new Benders cut
        if SolverResultsSbp.solver.termination_condition == TerminationCondition.infeasible:
            # the problem has to be feasible because I am not able to obtain the sum of infeasibilities of the phase I
            Delta[l] =  0
        else:
            # updating lower and upper bound
            Z_Lower =              Z1
            Z_Upper = min(Z_Upper, Z1 - pyo.value(mMaster_Bd.vTheta) + Z2)
            print('Iteration ', l, ' Z_Lower ... ', Z_Lower)
            print('Iteration ', l, ' Z_Upper ... ', Z_Upper)

            mMaster_Bd.vTheta.free()

            Delta[l] =  1

        for i,j in mFCTP.i*mFCTP.j:
            PI_L[l,i,j] = mFCTP.dual[mFCTP.eFlowLimit[i,j]]

        mMaster_Bd.vY.unfix()

        # add one cut
        mMaster_Bd.ll.add(l)
        ll = mMaster_Bd.ll
        mMaster_Bd.eBd_Cuts = Constraint(mMaster_Bd.ll, rule=eBd_Cuts, doc='Benders cuts')

# mFCTP.eCostSubp.deactivate()
# mFCTP.vY.unfix()
#
# def eCost(mFCTP):
#     return sum(mFCTP.pF[i,j]*mFCTP.vY[i,j] for i,j in mFCTP.i*mFCTP.j) + sum(mFCTP.pC[i,j]*mFCTP.vX[i,j] for i,j in mFCTP.i*mFCTP.j) + sum(mFCTP.vDNS[j]*1000 for j in mFCTP.j)
# mFCTP.eCost = Objective(rule=eCost, sense=minimize, doc='total cost')
#
# SolverResults = Solver.solve(mFCTP, tee=True)
# SolverResults.write()
