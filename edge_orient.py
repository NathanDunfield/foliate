import find_orient, link, dual_cellulation, surface
import snappy
import snappy.snap.t3mlite as t3m
from snappy.snap.t3mlite.simplex import (Head, Tail,
                                         ZeroSubsimplices, OneSubsimplices,
                                         RightFace, LeftFace)
from sage.all import Graph

class EdgeOrientation(object):
    """
    An orientation on the edges of a triangulation of a closed
    3-manifold where no face is a directed cycle.

    >>> N = t3m.Mcomplex('jLLvQPQcdfhghigiihshhgfifme')
    >>> orients = edge_orientations(N)
    >>> good = [eo for eo in orients if eo.link_compatible_with_foliation()]
    >>> len(good)
    8
    >>> [eo(e) for e in N.Edges] == eo.signs
    True
    >>> tet = N.Tetrahedra[3]
    >>> [eo.is_long(tet, e) for e in OneSubsimplices]
    [False, False, True, False, False, False]
    >>> eo.has_super_long_edge()
    False
    >>> eo.gives_foliation()
    True
    """
    def __init__(self, mcomplex, link_sphere, signs):
        self.mcomplex, self.signs = mcomplex, signs
        self.vertex_link = link_sphere
        assert len(mcomplex.Vertices) == 1 and mcomplex.Vertices[0].link_genus() == 0
        self._add_link_vertex_signs()

    def __call__(self, edge):
        """
        Return the sign self assigns to the given edge.
        """
        return self.signs[edge.Index]

    def _add_link_vertex_signs(self):
        self.link_vertex_signs = dict()
        for vert in self.vertex_link.vertices:
            i = vert.index
            edge_sign = self.signs[abs(i) - 1]
            vert_sign = 1 if edge_sign*i > 0 else -1
            self.link_vertex_signs[vert] = vert_sign

    def link_subgraphs(self):
        pos, neg = [], []
        for vert in self.vertex_link.vertices:
            if self.link_vertex_signs[vert] > 0:
                pos.append(vert.index)
            else:
                neg.append(vert.index)
        assert {abs(p) for p in pos} == {abs(n) for n in neg}
        G = self.vertex_link.edge_graph()
        return G.subgraph(vertices=pos), G.subgraph(vertices=neg)

    def link_compatible_with_foliation(self):
        pos, neg = self.link_subgraphs()
        return pos.is_connected() and neg.is_connected()

    def local_structure(self, tet, vertex):
        """
        Returns the number of "out" and "in" arrows from the given vertex.
        """
        a = vertex
        signs = []
        for b in ZeroSubsimplices:
            if b != a:
                e = tet.Class[a|b]
                signs.append(self(e) * e.orientation_with_respect_to(tet, a, b))
        return signs.count(1), signs.count(-1)

    def local_structure_edge(self, tet, edge):
        return sorted([self.local_structure(tet, Head[edge]),
                       self.local_structure(tet, Tail[edge])])

    def is_long(self, tet_or_corner, edge=None):
        if edge is None:  # Called with a corner
            tet = tet_or_corner.Tetrahedron
            edge = tet_or_corner.Subsimplex
        else:
            tet = tet_or_corner

        for a in [Head[edge], Tail[edge]]:
            if 0 not in self.local_structure(tet, a):
                return False
        return True
                        
    def has_super_long_edge(self):
        """
        A super long edge is an edge which is the long edge in every
        adjacent tetrahedra.  By definition, the long edge in a
        tetrahedron is the unique one the runs from the source to the
        sink.
        """
        for edge in self.mcomplex.Edges:
            if all(self.is_long(c) for c in edge.Corners):
                return True
        return False

    def num_sutures(self):
        M = self.mcomplex
        G = Graph()
        for tet in M.Tetrahedra:
            for e in OneSubsimplices:
                if self.local_structure_edge(tet, e) in [[(0,3), (1,2)], [(2,1), (3,0)]]:
                    i = tet.Class[RightFace[e]].Index
                    j = tet.Class[LeftFace[e]].Index
                    G.add_edge( (i,j) )
        assert G.num_verts() == len(M.Faces)
        return len(G.connected_components())

    def gives_foliation(self):
        if self.has_super_long_edge():
            return False

        ans1 = self.link_compatible_with_foliation()
        ans2 = self.num_sutures() == 1
        assert ans1 == ans2
        return ans1
    

class IdealEdgeOrientation(EdgeOrientation):
    """
    An orientation on the edges of an ideal triangulation of a
    1-cusped 3-manifold where no face is a directed cycle.

    >>> N = t3m.Mcomplex('v1234')
    >>> orients = edge_orientations(N)
    >>> [eo.num_sutures() for eo in orients]
    [2, 2]
    """
    def __init__(self, mcomplex, link_triangulation, link_dual_cellulation, signs):
        self.mcomplex, self.signs = mcomplex, signs
        self.vertex_link = link_triangulation
        self.link_dual_cellulation = link_dual_cellulation
        assert len(mcomplex.Vertices) == 1 and mcomplex.Vertices[0].link_genus() == 1
        self._add_link_vertex_signs()
        
    def sutures(self):
        """
        Returns a list of cycles on the 1-skeleton of the dual cellulation
        of the vertex link, each representing one of the sutures
        induced by the corresponding branched surface.

        >>> N = t3m.Mcomplex('m015')
        >>> orients = edge_orientations(N)
        >>> len([eo.sutures() for eo in orients])
        2
        """
        T = self.vertex_link
        D = self.link_dual_cellulation
        weights = len(T.edges) * [0]
        for tri in T.triangles:
            signs = [self.link_vertex_signs[v] for v in tri.vertices]
            if len(set(signs)) > 1:
                if sorted(signs) == [-1, 1, 1]:
                    i = signs.index(-1)
                    j = (i + 1) % 3
                    side = surface.Side(tri, (i, j))
                else:
                    assert sorted(signs) == [-1, -1, 1]
                    i = signs.index(1)
                    j = (i - 1) % 3
                    side = surface.Side(tri, (j, i))
                edge = side.edge()
                sign = edge.orientation_with_respect_to(side)
                # Edge and its dual have the same index
                weights[edge.index] = sign
        return dual_cellulation.OneCycle(D, weights).components()
                

    def link_compatible_with_foliation(self):
        pass

def edge_orientations(manifold):
    if isinstance(manifold, t3m.Mcomplex):
        N = manifold
        M = None
    else: # SnapPy manifold
        M = manifold
        N = t3m.Mcomplex(M)        
    if len(N.Vertices) == 1 and N.Vertices[0].link_genus() == 0:
        vertex_link = link.LinkSphere(N)
        for signs in find_orient.cycle_free_orientations(N):
            yield EdgeOrientation(N, vertex_link, signs)
    else:
        if M is not None:
            N, vertex_link, dual_cell, (mstar, lstar) = dual_cellulation.peripheral_curve_package(M)
            dual_cell.meridian_star = mstar
            dual_cell.longitude_star = lstar
        else:
            vertex_link = link.LinkSurface(N)
            dual_cell  = dual_cellulation.DualCellulation(vertex_link)
        assert len(N.Vertices) == 1 and N.Vertices[0].link_genus() == 1
        for signs in find_orient.cycle_free_orientations(N):
            yield IdealEdgeOrientation(N, vertex_link, dual_cell, signs)

def test_sutures(n=100, progress=True):
    """
    >>> test_sutures(5, False)
    """
    census = snappy.OrientableCuspedCensus(cusps=1)
    for i in range(n):
        M = census.random()
        sutures = [eo.sutures() for eo in edge_orientations(M)]
        if len(sutures):
            D = sutures[0][0].cellulation
            mstar, lstar = D.meridian_star, D.longitude_star
            sutures = [[(mstar(s), lstar(s)) for s in suture] for suture in sutures]
        if progress:
            print(M.name() + ' ' + repr(sutures))

        
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    #N = t3m.Mcomplex('jLLvQPQcdfhghigiihshhgfifme')
    #orients = edge_orientations(N)
    #[eo.num_sutures() for eo in orients]
    N = t3m.Mcomplex('m004')
    orients = edge_orientations(N)
    eo = next(orients)
    C = eo.vertex_link
    s = eo.sutures()
