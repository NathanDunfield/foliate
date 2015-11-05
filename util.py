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
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()

        
