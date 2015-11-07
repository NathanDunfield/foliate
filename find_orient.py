"""
Given a triangulation of a 3-manifold, find all orientations of
the one-skeleton where no triangular face is a directed cycle.
"""

import snappy
import snappy.snap.t3mlite as t3m
from snappy.snap.t3mlite.simplex import *
from sage.sat.solvers import CryptoMiniSat
from sage.all import infinity

# -------- t3m preliminaries --------

# A simplex is oriented like this:  
#     1     
#    /|\    
#   / | \   
#  /  |  \  
# 2---|---3 
#  \  |  /  
#   \ | /   
#    \|/    
#     0
#
#
# where the outward normal induces the following anticlockwise
# orientation on the 2-faces.

VerticesOfFace = { F0 : (V1, V3, V2), F1 : (V0, V2, V3),
                   F2 : (V0, V3, V1), F3 : (V0, V1, V2) }

def oriented_edges_around_faces(triangulation):
    ans = []
    for tet in triangulation.Tetrahedra:        
        for vertices in VerticesOfFace.values():
            face = []
            for i in range(3):
                a, b = vertices[i], vertices[(i+1)%3]
                edge = tet.Class[a|b]
                sign = edge.orientation_with_respect_to(tet, a, b)
                face.append(sign*(edge.Index + 1))
            ans.append(face)
    return ans

def all_solutions(solver):
    """
    Return all solutions of a CryptoMiniSat solver. 

    Note: Modifies the solver inplace.
    """
    while True:
        solution = solver()
        if solution == False:
            return
        else:
            yield solution
            # Add a clause which excludes the solution just found.
            clause = [-i if s else i for i, s in enumerate(solution)]
            solver.add_clause(tuple(clause[1:]))

def cycle_free_orientations(triangulation):
    """
    Returns all orientations of the one-skeleton where no triangular
    face is a directed cycle.  Orientations are given relative to the
    default orientation as a sequence of 1's and -1's.

    >>> M = t3m.Mcomplex('jLvMLQQbfefgihhiixiptvvvgof')
    >>> list(cycle_free_orientations(M))
    []
    >>> M = t3m.Mcomplex('jLvLQAQbffghghiiieuaiikktuu')
    >>> len(list(cycle_free_orientations(M)))
    10
    """
    solver = CryptoMiniSat()
    for face_data in oriented_edges_around_faces(triangulation):
        solver.add_clause(face_data)
    # By symmetry, we might as well assume the first edge is
    # positively oriented
    solver.add_clause((1,))
    for sol in all_solutions(solver):
        yield [1 if s else -1 for s in sol[1:]]
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()

        
