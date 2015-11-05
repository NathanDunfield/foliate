import edge_orient
import snappy
import snappy.snap.t3mlite as t3m
import pandas as pd

def has_compatible_foliation(snappy_manifold):
    ans = first_foliation(snappy_manifold)
    return ans is not None

def first_foliation(snappy_manifold):
    T = t3m.Mcomplex(snappy_manifold.filled_triangulation())
    if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
        for eo in edge_orient.edge_orientations(T):
            if eo.gives_foliation():
                return eo

def make_file():
    file = open('/tmp/foliation_info.csv', 'w')
    file.write('name,known_fol\n')
    for M in snappy.OrientableClosedCensus:
        file.write('"%s",%s\n' % (M, has_compatible_foliation(M)))
        file.flush()

def disorder(snappy_manifold):
    pass

def compare_data():
    dh = pd.read_csv('/Users/dunfield/h/derived_data/combined_data_2014_3_12.csv')
    df = pd.read_csv('/tmp/foliation_info.csv')
    da = dh.merge(df, on='name')[['name', 'L_space', 'known_fol']]
    bad = da[(da.L_space==1)&(da.known_fol)]
    return da, bad

def repeatibility():
    for M in snappy.OrientableClosedCensus[:100]:
        print M, {M.filled_triangulation().triangulation_isosig() for i in range(100)}
            
