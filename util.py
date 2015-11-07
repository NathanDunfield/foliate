import snappy
import snappy.snap.t3mlite as t3m


# Function not used, but kept around for doctests
def closed_from_isosig(isosig):
    """
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
                ans.add(T.triangulation_isosig())
            N.randomize()

    return sorted(ans, key=lambda i:i.swapcase())
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
