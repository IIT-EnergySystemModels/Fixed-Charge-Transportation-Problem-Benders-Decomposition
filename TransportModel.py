# Developed by
#
#    Andres Ramos
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingenieria - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu
#    https://pascua.iit.comillas.edu/aramos/Ramos_CV.htm
#
#    May 26, 2025

import pyomo.environ as pyo
from   pyomo.environ import ConcreteModel, Set, Param, Var, NonNegativeReals, Constraint, Objective, minimize, Suffix
from   pyomo.opt     import SolverFactory

mTransport = ConcreteModel('Transportation Problem')

mTransport.i = Set(initialize=['Vigo',   'Algeciras'            ], doc='origins'     )
mTransport.j = Set(initialize=['Madrid', 'Barcelona', 'Valencia'], doc='destinations')

mTransport.pA = Param(mTransport.i, initialize={'Vigo'  : 350, 'Algeciras': 700                 }, doc='origin capacity'   )
mTransport.pB = Param(mTransport.j, initialize={'Madrid': 400, 'Barcelona': 450, 'Valencia': 150}, doc='destination demand')

TransportationCost = {
    ('Vigo',      'Madrid'   ): 0.06,
    ('Vigo',      'Barcelona'): 0.12,
    ('Vigo',      'Valencia' ): 0.09,
    ('Algeciras', 'Madrid'   ): 0.05,
    ('Algeciras', 'Barcelona'): 0.15,
    ('Algeciras', 'Valencia' ): 0.11,
    }

mTransport.pC = Param(mTransport.i*mTransport.j, initialize=TransportationCost, doc='per unit transportation cost')

mTransport.vX = Var  (mTransport.i*mTransport.j, bounds=(0.0,None), doc='units transported', within=NonNegativeReals)

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

Solver = SolverFactory('gurobi')
Solver.options['LogFile'] = 'mTransport.log'
SolverResults = Solver.solve(mTransport, tee=True)
SolverResults.write()
mTransport.pprint()

mTransport.vX.display()
for j in mTransport.j:
    print(mTransport.dual[mTransport.eDemand[j]])
