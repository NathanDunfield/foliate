import edge_orient
import snappy
import snappy.snap.t3mlite as t3m
import pandas as pd
import util

def has_compatible_foliation(snappy_manifold):
    ans = first_foliation(snappy_manifold)
    return ans is not None

def first_foliation_mcomplex(mcomplex):
    for eo in edge_orient.edge_orientations(mcomplex):
        if eo.gives_foliation():
            return eo
    
def first_foliation(snappy_manifold, rand_max, examine_max):
    for iso in util.closed_isosigs(snappy_manifold, rand_max)[:examine_max]:
        T = t3m.Mcomplex(iso)
        T.name = iso
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orient = list(edge_orient.edge_orientations(T))
            for eo in orient:
                if eo.gives_foliation():
                    return eo

def nonorderable(snappy_manifold, max_triangulations=10):
    for iso in util.closed_isosigs(snappy_manifold)[:max_triangulations]:
        print iso
        T = t3m.Mcomplex(iso)
        T.name = iso
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orients = edge_orient.edge_orientations(T)
            found = False
            for eo in orients:
                found = True
                break
            if not found:
                return iso
                    
            

def make_fol_file():
    file = open('/tmp/foliation_info.csv', 'w')
    file.write('name,taut,good_tri\n')
    for M in snappy.OrientableClosedCensus:
        fol = first_foliation(M, 10)
        if fol is None:
            file.write('"%s",0,\n' % M)
        else:
            file.write('"%s",1,"%s"\n' % (M, fol.mcomplex.name))
        file.flush()

def make_disorder_file():
    file = open('/tmp/orderable_info.csv', 'w')
    file.write('name,orderable,good_tri\n')
    for M in snappy.OrientableClosedCensus:
        ans = nonorderable(M, 100)
        if ans is None:
            file.write('"%s",0,\n' % M)
        else:
            file.write('"%s",-1,"%s"\n' % (M, ans))
        file.flush()


def compare_data():
    dh = pd.read_csv('/Users/dunfield/h/derived_data/combined_data_2014_3_12.csv')
    df = pd.read_csv('foliation_info.csv')
    da = dh.merge(df, on='name')[['name', 'L_space', 'taut']]
    bad = da[(da.L_space==1)&(da.taut)]
    assert len(bad) == 0
    return da

def repeatibility():
    for M in snappy.OrientableClosedCensus[:100]:
        print M, {M.filled_triangulation().triangulation_isosig() for i in range(100)}
            
