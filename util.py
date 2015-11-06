import snappy, regina
import snappy.snap.t3mlite as t3m

def closed_from_file(filename):
    data = t3m.files.read_SnapPea_file(filename)
    tets = t3m.mcomplex.tets_from_data(data)
    return t3m.Mcomplex(tets)

def closed_from_isosig(isosig):
    """
    >>> isosig = 'jLLvQPQcdfhghigiihshhgfifme'
    >>> N = closed_from_isosig(isosig)
    >>> len(N)
    9
    >>> N.snappy_manifold().identify()
    [m003(-2,3)]
    """
    T = regina.NTriangulation(isosig)
    data = t3m.files.read_SnapPea_file(data=T.snapPea())
    tets = t3m.mcomplex.tets_from_data(data)
    return t3m.Mcomplex(tets)
        
def closed_isosigs(snappy_manifold, trys=20):
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
            ans.add(T.triangulation_isosig())
            N.randomize()

    return sorted(ans, key=lambda i:i.swapcase())
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
