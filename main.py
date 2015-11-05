import orient, link
import snappy
import snappy.snap.t3mlite as t3m

def has_compatible_foliation(snappy_manifold):
    ans = first_foliation(snappy_manifold)
    return ans is not None

def first_foliation(snappy_manifold):
    T = t3m.Mcomplex(snappy_manifold.filled_triangulation())
    if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
        S = link.LinkSphere(T)
        for prefol in orient.cycle_free_orientations(T):
            print snappy_manifold, S.gives_foliation(prefol)
            if S.gives_foliation(prefol):
                return T, S, prefol

def make_file():
    file = open('/tmp/foliation_info.csv', 'w')
    file.write('name,known_fol\n')
    for M in snappy.OrientableClosedCensus[:100]:
        file.write('"%s",%s\n' % (M, has_compatible_foliation(M)))
        file.flush()

def disorder(snappy_manifold):
    pass
    

def repeatibility():
    for M in snappy.OrientableClosedCensus[:100]:
        print M, {M.filled_triangulation().triangulation_isosig() for i in range(100)}
            
