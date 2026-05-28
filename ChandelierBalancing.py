# Developed by
#    Andres Ramos
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingenieria - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu

# May 27, 2025

import pandas as pd
from pyomo.environ import ConcreteModel, Binary, Set, Param, Var, Constraint
from pyomo.opt     import SolverFactory

SolverName = 'gurobi'

# model declaration
mChB = ConcreteModel('Chandelier Balancing')


def plain_run(mChB, SolverName):

    mChB.p = Set(initialize=['a',  'b',  'c',  'd',  'e',  'f',  'g',  'h',  'i' ], ordered=True, doc='positions')
    mChB.w = Set(initialize=['w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9'], ordered=True, doc='weights'  )

    # parameters
    mChB.weights = Param(mChB.w, initialize={'w1': 1, 'w2': 2, 'w3': 3, 'w4': 4, 'w5': 5, 'w6': 6, 'w7': 7, 'w8': 8, 'w9': 9}, doc='weight of each mass')

    # Variables
    mChB.location = Var(mChB.w*mChB.p, within=Binary, doc='locate weight w in position p')

    def branchABC(mChB):
        return sum(mChB.location[w,'a']*mChB.weights[w] for w in mChB.w)*2 - sum(mChB.location[w,'b']*mChB.weights[w] for w in mChB.w) - sum(mChB.location[w,'c']*mChB.weights[w] for w in mChB.w)*2 == 0
    mChB.branchABC = Constraint(rule=branchABC, doc='balance of branch ABC')

    def branchDEF(mChB):
        return sum(mChB.location[w,'d']*mChB.weights[w] for w in mChB.w)*2 + sum(mChB.location[w,'e']*mChB.weights[w] for w in mChB.w) - sum(mChB.location[w,'f']*mChB.weights[w] for w in mChB.w)   == 0
    mChB.branchDEF = Constraint(rule=branchDEF, doc='balance of branch DEF')

    def branchGHI(mChB):
        return sum(mChB.location[w,'g']*mChB.weights[w] for w in mChB.w)*2 + sum(mChB.location[w,'h']*mChB.weights[w] for w in mChB.w) - sum(mChB.location[w,'i']*mChB.weights[w] for w in mChB.w)*3 == 0
    mChB.branchGHI = Constraint(rule=branchGHI, doc='balance of branch GHI')

    def branchTot(mChB):
        return ((sum(mChB.location[w,'a']*mChB.weights[w] for w in mChB.w) + sum(mChB.location[w,'b']*mChB.weights[w] for w in mChB.w) + sum(mChB.location[w,'c']*mChB.weights[w] for w in mChB.w))*3 -
                (sum(mChB.location[w,'d']*mChB.weights[w] for w in mChB.w) + sum(mChB.location[w,'e']*mChB.weights[w] for w in mChB.w) + sum(mChB.location[w,'f']*mChB.weights[w] for w in mChB.w))*2 -
                (sum(mChB.location[w,'g']*mChB.weights[w] for w in mChB.w) + sum(mChB.location[w,'h']*mChB.weights[w] for w in mChB.w) + sum(mChB.location[w,'i']*mChB.weights[w] for w in mChB.w))*3) == 0
    mChB.branchTot = Constraint(rule=branchTot, doc='total balance')

    def obligation(mChB, w):
        return sum(mChB.location[w,p] for p in mChB.p) == 1
    mChB.obligation = Constraint(mChB.w, rule=obligation, doc='weights must be located in one position')

    def oneweightp(mChB, p):
        return sum(mChB.location[w,p] for w in mChB.w) == 1
    mChB.oneweightp = Constraint(mChB.p, rule=oneweightp, doc='one weight by position')

    Solver = SolverFactory(SolverName)
    mChB.write('chandelier.lp', io_options={'symbolic_solver_labels': True})
    SolverResults = Solver.solve(mChB, tee=True)
    SolverResults.write()  # summary of the solver results
    mChB.location.pprint()

if __name__ == "__main__":
    plain_run(mChB, SolverName)
