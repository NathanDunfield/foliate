"""
The dual cellulation to a triangulation of an oriented surface.  
"""

from snappy.snap.t3mlite import Mcomplex
from link import LinkSurface
import surface
import snappy
from sage.all import (ZZ, matrix, vector, ChainComplex, cached_method)

class DualCell(object):
    """
    A cell in the dual cellulation
    """
    def __init__(self, dual_cell):
        self.dual_cell = dual_cell
        self.index = dual_cell.index

class Vertex(DualCell):
    """
    A vertex of the dual cellulation.
    """        

class Edge(DualCell):
    """
    An oriented edge of the dual cellulation. The orientation
    convention is that if e is the original edge and d is the dual
    edge then one rotates from d to e anticlockwise.
    """
    def __init__(self, dual_cell):
        DualCell.__init__(self, dual_cell)
        self.vertices = [None, None]
    
class Face(DualCell):
    """
    A face of the dual cellulation, which is just an n-gon.  The
    vertices are numbered range(n) in an anti-clockwise orientation.

    Oriented edges of the face are specified by their initial vertex.
    """
    def __init__(self, dual_cell):
        DualCell.__init__(self, dual_cell)
        self.n = len(dual_cell.corners)
        self.edges_with_orientations = self.n*[None]

    def __repr__(self):
        return "<Face: %s>" % self.index

    #def oriented_sides(self):
    #    return [Side(self, e) for e in oriented_edges_of_triangle]

class DualCellulation(object):
    """
    The dual cellulation to a triangulation of a surface. 
    """
    def __init__(self, triangulation):
        self.dual_triangulation = triangulation
        self.from_original = dict()
        self.vertices = []
        self.edges = []
        self.faces = []
        for vertex in triangulation.vertices:
            face = Face(vertex)
            self.faces.append(face)
            self.from_original[vertex] = face
        for edge in triangulation.edges:
            dual_edge = Edge(edge)
            self.edges.append(dual_edge)
            self.from_original[edge] = dual_edge
        for triangle in triangulation.triangles:
            vertex = Vertex(triangle)
            self.vertices.append(vertex)
            self.from_original[triangle] = vertex

        for vertex in triangulation.vertices:
            face = self.from_original[vertex]
            for i, corner in enumerate(vertex.corners):
                edge, orient = corner.edge_with_orientation()
                dual_edge = self.from_original[edge]
                face.edges_with_orientations[i] = (dual_edge, -orient)
                if orient > 0:
                    dual_edge.vertices = [self.from_original[side.triangle] for side in edge.sides]        
                
    def euler(self):
        """
        >>> N = Mcomplex('o9_12345')
        >>> D = DualCellulation(LinkSurface(N))
        >>> D.euler()
        0
        
        >>> N = Mcomplex('jLLvQPQcdfhghigiihshhgfifme')
        >>> D = DualCellulation(LinkSurface(N))
        >>> D.euler()
        2
        """
        return len(self.vertices) - len(self.edges) + len(self.faces)

    @cached_method
    def B1(self):
        """
        The matrix describing the boundary map C_1 -> C_0
        """
        V, E = len(self.vertices), len(self.edges)
        assert range(V) == sorted(v.index for v in self.vertices)
        assert range(E) == sorted(e.index for e in self.edges)
        
        D = matrix(ZZ, V, E, sparse=True)
        for e in self.edges:
            v_init = e.vertices[0].index
            v_term = e.vertices[1].index
            D[v_term, e.index] += 1
            D[v_init, e.index] += -1

        return D

    @cached_method
    def B2(self):
        """
        The matrix describing the boundary map C_2 -> C_1

        Does *not* assume that the faces are numbered like
        range(len(faces)).
        """
        E, F = len(self.edges), len(self.faces)
        assert range(E) == sorted(e.index for e in self.edges)
        D = matrix(ZZ, E, F, sparse=True)
        for i, face in enumerate(self.faces):
            for edge, sign in face.edges_with_orientations:
                D[edge.index, i] += sign

        return D

    @cached_method
    def chain_complex(self):
         return ChainComplex( {1:self.B1(), 2:self.B2()} , degree=-1 )

    def homology_test(self):
        T = self.dual_triangulation
        B1, B2 = self.B1(), self.B2()
        assert B1*B2 == 0
        assert T.euler() == self.euler()
        CD = self.chain_complex()
        CT = T.chain_complex()
        assert CD.homology() == CT.homology()

def test_homologically(n=100, progress=True):
    """
    >>> test_homologically(5, False)
    """
    for i in range(n):
        M = snappy.HTLinkExteriors.random()
        N = Mcomplex(M)
        C = LinkSurface(N)
        C.index()   # So that homology works
        D = DualCellulation(C)
        D.homology_test()
        if progress:
            H = D.chain_complex().homology()
            print(M.name() + ' ' + repr(H))

                

          
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    N = Mcomplex('s000')
    C = LinkSurface(N)
    C.index()
    D = DualCellulation(C)
    D.homology_test()
