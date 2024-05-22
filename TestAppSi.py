import pyomo.environ as pyo
from   pyomo.environ import ConcreteModel, Set, Param, Var, NonNegativeReals, Constraint, Objective, minimize, Suffix
from   pyomo.opt     import SolverFactory
from pyomo.contrib import appsi
import numpy as np
from pyomo.common.timing import HierarchicalTimer

mTransport = ConcreteModel('Transportation Problem')

mTransport.i = Set(initialize=['Vigo',   'Algeciras'            ], doc='origins'     )
mTransport.j = Set(initialize=['Madrid', 'Barcelona', 'Valencia'], doc='destinations')

mTransport.pA = Param(mTransport.i, initialize={'Vigo'  : 350, 'Algeciras': 700                 }, doc='origin capacity'   )
mTransport.pB = Param(mTransport.j, initialize={'Madrid': 400, 'Barcelona': 450, 'Valencia': 150}, doc='destination demand', mutable=True)
mTransport.pD = Param(mTransport.j, initialize={'Madrid': 400, 'Barcelona': 450, 'Valencia': 150}, doc='destination demand', mutable=True)

TransportationCost = {
    ('Vigo',      'Madrid'   ): 0.06,
    ('Vigo',      'Barcelona'): 0.12,
    ('Vigo',      'Valencia' ): 0.09,
    ('Algeciras', 'Madrid'   ): 0.05,
    ('Algeciras', 'Barcelona'): 0.15,
    ('Algeciras', 'Valencia' ): 0.11,
    }

mTransport.pC = Param(mTransport.i, mTransport.j, initialize=TransportationCost, doc='per unit transportation cost')

mTransport.vX = Var  (mTransport.i, mTransport.j, bounds=(0.0,None), doc='units transported', within=NonNegativeReals)

def eCapacity(mTransport, i):
    return sum(mTransport.vX[i,j] for j in mTransport.j) <= mTransport.pA[i]
mTransport.eCapacity = Constraint(mTransport.i, rule=eCapacity, doc='maximum capacity of each origin')

def eDemand  (mTransport, j):
    return sum(mTransport.vX[i,j] for i in mTransport.i) >= mTransport.pB[j]
mTransport.eDemand   = Constraint(mTransport.j, rule=eDemand,   doc='demand supply at destination'   )

def eCost(mTransport):
    return sum(mTransport.pC[i,j]*mTransport.vX[i,j] for i,j in mTransport.i*mTransport.j)
mTransport.eCost = Objective(rule=eCost, sense=minimize, doc='transportation cost')

mTransport.write('mTransport.lp', io_options={'symbolic_solver_labels': True})

mTransport.dual = Suffix(direction=Suffix.IMPORT)

opt = appsi.solvers.Gurobi()
timer = HierarchicalTimer()
for p_val in np.linspace(0.8, 1, num=3):
    for j in mTransport.j:
        mTransport.pB[j] = float(p_val)*mTransport.pD[j]
    res = opt.solve(mTransport, timer=timer)
    assert res.termination_condition == appsi.base.TerminationCondition.optimal
    print(mTransport.eCost.expr())

    mTransport.vX.display()
    mTransport.eDemand.display()
    # for j in mTransport.j:
    print(opt.get_duals(mTransport.eDemand[:]).values())

# print(timer)