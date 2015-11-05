import find_orient, link
from snappy.snap.t3mlite.simplex import Head, Tail, ZeroSubsimplices, OneSubsimplices
from util import closed_from_isosig

class EdgeOrientation(object):
    """
    An orientation on the edges of a triangulation of a closed
    3-manifold where no face is a directed cycle.

    >>> N = closed_from_isosig('jLLvQPQcdfhghigiihshhgfifme')
    >>> orients = edge_orientations(N)
    >>> good = [eo for eo in orients if eo.link_compatible_with_foliation()]
    >>> len(good)
    8
    >>> [eo(e) for e in N.Edges] == eo.signs
    True
    >>> tet = N.Tetrahedra[3]
    >>> [eo.is_long(tet, e) for e in OneSubsimplices]
    [False, False, False, False, True, False]
    >>> eo.has_super_long_edge()
    False
    >>> eo.gives_foliation()
    True
    """
    def __init__(self, mcomplex, link_sphere, signs):
        self.mcomplex, self.signs = mcomplex, signs
        self.link_sphere = link_sphere

    def __call__(self, edge):
        """
        Return the sign self assigns to the given edge.
        """
        return self.signs[edge.Index]

    def link_subgraphs(self):
        pos, neg = [], []
        for vert in self.link_sphere.vertices:
            i = vert.index
            s = self.signs[abs(i) - 1]
            if s*i > 0:
                pos.append(i)
            else:
                neg.append(i)
        assert {abs(p) for p in pos} == {abs(n) for n in neg}
        G = self.link_sphere.edge_graph()
        return G.subgraph(vertices=pos), G.subgraph(vertices=neg)

    def link_compatible_with_foliation(self):
        pos, neg = self.link_subgraphs()
        return pos.is_connected() and neg.is_connected()

    def is_long(self, tet_or_corner, edge=None):
        if edge is None:  # Called with a corner
            tet = tet_or_corner.Tetrahedron
            edge = tet_or_corner.Subsimplex
        else:
            tet = tet_or_corner

        for a in [Head[edge], Tail[edge]]:
            adjacent_signs = set()
            for b in ZeroSubsimplices:
                if b != a:
                    e = tet.Class[a|b]
                    adjacent_signs.add(self(e) * e.orientation_with_respect_to(tet, a, b))
            if len(adjacent_signs) != 1:
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

    def gives_foliation(self):
        """
        WARNING: Not yet proved this works!
        """
        return self.link_compatible_with_foliation() and not self.has_super_long_edge()
    
def edge_orientations(mcomplex):
    link_sphere = link.LinkSphere(mcomplex)
    for signs in find_orient.cycle_free_orientations(mcomplex):
        yield EdgeOrientation(mcomplex, link_sphere, signs)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    N = closed_from_isosig('jLLvQPQcdfhghigiihshhgfifme')
