import snappy
import snappy.snap.t3mlite as t3m
import pyregina, regina

def to_iso(snappy_manifold):
    M = snappy_manifold
    return (M.num_tetrahedra(), M.triangulation_isosig(decorated=False))
        
def fancy_closed_isosigs(snappy_manifold, height):
    """
    >> M = snappy.Manifold('m004(1,2)')
    >> len(fancy_closed_isosigs(M,0)) > 0
    True
    """
    M = snappy_manifold.copy()
    assert M.cusp_info('complete?') == [False]
    starts = set()

    for i in range(10):
        starts.add(to_iso(M.filled_triangulation()))
        for curve in M.dual_curves():
            N = M.drill(curve)
            N.dehn_fill((1,0), 1)
            starts.add(to_iso(N.filled_triangulation()))
        M.randomize()

    min_tets = min(starts)[0]
    max_tets = min_tets + height
    starts = sorted(N for N in starts if N[0] <= max_tets)
    isosigs = set()
    
    for n, iso in starts:
        if (n, iso) not in isosigs:
            print('Starting massive search from %s' % iso)
            T = pyregina.Triangulation(iso)
            isosigs.update(T.retriangulate(max_tets - n))
        else:
            print('Already found %s' % iso)
    return sorted(isosigs)

def two_vertex_tris(snappy_manifold, height=0):
    M = snappy_manifold.filled_triangulation()
    T = regina.NTriangulation(M._to_string())
    t = T.tetrahedra()[0]
    T.oneFourMove(t)
    isosig = T.isoSig()
    P = pyregina.Triangulation(isosig)
    isosigs = P.retriangulate(height)
    return isosigs

        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
