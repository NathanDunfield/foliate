import snappy
import snappy.snap.t3mlite as t3m


def closed_from_isosig(isosig):
    """
    Function not used, but kept around for doctests

    >>> isosig = 'jLLvQPQcdfhghigiihshhgfifme'
    >>> N = t3m.Mcomplex(isosig)
    >>> len(N)
    9
    >>> N.snappy_manifold().identify()
    [m003(-2,3)]
    >>> N.isosig() == isosig
    True
    """
    return t3m.Mcomplex(isosig)
        
def closed_isosigs(snappy_manifold, trys=20):
    """
    >>> M = snappy.Manifold('m004(1,2)')
    >>> len(closed_isosigs(M, trys=5)) > 0
    True
    """
    M = snappy_manifold.copy()
    assert M.cusp_info('complete?') == [False]
    surgery_descriptions = [M]

    for curve in M.dual_curves():
        N = M.drill(curve)
        N.dehn_fill((1,0), 1)
        surgery_descriptions.append(N.filled_triangulation([0]))

    ans = set()
    for N in surgery_descriptions:
        for i in range(trys):
            T = N.filled_triangulation()
            if T._num_fake_cusps() == 1:
                ans.add((T.num_tetrahedra(), T.triangulation_isosig(decorated=False)))
            N.randomize()

    return [iso for n, iso in sorted(ans)]

def cusped_triangulations(snappy_manifold, trys=1000):
    """
    >>> M = snappy.Manifold('m004')
    >>> len(cusped_triangulations(M, trys=100))
    1
    """
    M = snappy.Triangulation(snappy_manifold)
    ans = [M.copy()]
    for i in range(trys):
        M.randomize()
        if {len(M.isomorphisms_to(N)) > 0 for N in ans} == {False}:
            ans.append(M.copy())
    return ans
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
