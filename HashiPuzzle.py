# Developed by
#    Andres Ramos
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingenieria - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu

import pandas as pd
from pyomo.environ import ConcreteModel, NonNegativeIntegers, Set, Param, Var, Constraint, Objective, minimize
from pyomo.opt     import SolverFactory
from collections   import defaultdict

SolverName     = 'gurobi'

# model declaration
mHP = ConcreteModel('Hashi puzzle')


def plain_run(mHP, SolverName):

    mHP.x = Set(initialize=['x00', 'x01', 'x02', 'x03', 'x04', 'x05', 'x06', 'x07', 'x08', 'x09', 'x10', 'x11'], ordered=True, doc='abscises' )
    mHP.y = Set(initialize=['y00', 'y01', 'y02', 'y03', 'y04', 'y05', 'y06', 'y07', 'y08', 'y09', 'y10', 'y11'], ordered=True, doc='ordinates')

    Nodes = {
        ('x05', 'y01'): 1,
        ('x02', 'y02'): 3,
        ('x04', 'y02'): 3,
        ('x06', 'y02'): 4,
        ('x03', 'y03'): 7,
        ('x05', 'y03'): 6,
        ('x09', 'y03'): 1,
        ('x04', 'y04'): 7,
        ('x06', 'y04'): 6,
        ('x10', 'y04'): 7,
        ('x07', 'y05'): 5,
        ('x09', 'y05'): 7,
        ('x11', 'y05'): 4,
        ('x06', 'y06'): 5,
        ('x08', 'y06'): 6,
    }

    Arcs = set()
    Neighbors = defaultdict(list)
    for x,y,xx,yy in mHP.x*mHP.y*mHP.x*mHP.y:
        if (x,y) in Nodes and (xx,yy) in Nodes:
            if mHP.x.ord(x) > 2 and mHP.x.ord(x) < len(mHP.x)-2 and mHP.y.ord(y) > 1 and mHP.y.ord(y) < len(mHP.y)-1:
                if (xx,yy) == (mHP.x.next(x),mHP.y.next(y)) or (xx,yy) == (mHP.x.next(x,2),y) or (xx,yy) == (mHP.x.next(x),mHP.y.prev(y)) or (xx,yy) == (mHP.x.prev(x),mHP.y.prev(y)) or (xx,yy) == (mHP.x.prev(x,2),y) or (xx,yy) == (mHP.x.prev(x),mHP.y.next(y)):
                    # add the arc and its neighbors to the list
                    Arcs.add((x, y, xx, yy))
                    Arcs.add((xx, yy, x, y))
                    Neighbors[x,y].append((xx,yy))

    # parameters
    mHP.pNodes = Param(mHP.x, mHP.y, initialize=Nodes, doc='number of connections per node')

    # Variables
    mHP.vConnection = Var(Arcs, within=NonNegativeIntegers, doc='connections between two neighbor nodes')

    def eTotalConnections(mHP):
        return sum(mHP.vConnection[x,y,xx,yy] for x,y,xx,yy in Arcs)
    mHP.eTotalConnections = Objective(rule=eTotalConnections, sense=minimize, doc='total number of connections')

    def eConnections(mHP, x, y):
        if mHP.x.ord(x) > 2 and mHP.x.ord(x) < len(mHP.x)-2 and mHP.y.ord(y) > 1 and mHP.y.ord(y) < len(mHP.y)-1:
            return sum(mHP.vConnection[x,y,xx,yy] + mHP.vConnection[xx,yy,x,y] for xx,yy in Neighbors[x,y]) == mHP.pNodes[x,y]
        else:
            return Constraint.Skip
    mHP.eConnections = Constraint(Nodes, rule=eConnections, doc='connections on every node')

    Solver = SolverFactory(SolverName)
    mHP.write('connection.lp', io_options={'symbolic_solver_labels': True})
    SolverResults = Solver.solve(mHP, tee=True)
    SolverResults.write()  # summary of the solver results
    mHP.vConnection.pprint()

if __name__ == "__main__":
    plain_run(mHP, SolverName)
