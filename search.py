import edge_orient, find_orient
import snappy
import snappy.snap.t3mlite as t3m
import util
from collections import Counter
import regina
import random

def edge_orientation_stats(manifold):
    if isinstance(manifold, str):
        T = t3m.Mcomplex(manifold)
    else:
        T = manifold
    ans = dict()
    if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
        orient = list(edge_orient.edge_orientations(T))
        ans['num_orient'] = len(orient)
        ans['sutures'] = Counter(eo.num_sutures() for eo in orient)
        ans['num_sinks'] = Counter(eo.num_super_long_edges() for eo in orient)
        ans['fol'] = Counter(eo.gives_foliation() for eo in orient)
    return ans


def first_foliation(snappy_manifold, max_triangulations=10):
    for iso in util.closed_isosigs(snappy_manifold)[:max_triangulations]:
        T = t3m.Mcomplex(iso)
        T.name = iso
        if len(T.Vertices) == 1 and T.Vertices[0].link_genus() == 0:
            orient = list(edge_orient.edge_orientations(T))
            for eo in orient:
                if eo.gives_foliation():
                    return eo

def try_persistent(snappy_manifold):
    Y = snappy_manifold
    n = len(Y.dual_curves())
    for i in range(n):
        M = Y.drill(i).filled_triangulation()
        if M.solution_type() == 'all tetrahedra positively oriented': 
            slopes = edge_orient.degeneracy_slopes_with_search(M)
            print slopes

def tri_supports_foliation(iso):
    T = t3m.Mcomplex(iso)
    for eo in edge_orient.edge_orientations(T):
        if eo.gives_foliation():
            return eo
        

        
    
            
def highest_degree_vertex(regina_tri):
    return max(regina_tri.vertices(), key=lambda v:v.degree())

def triangles_on_high(regina_tri):
    T = regina_tri
    v = highest_degree_vertex(T)
    for F in T.triangles():
        if F.vertex(0) == F.vertex(1) == F.vertex(2) == v:
            return F

def obvious_simplify(regina_tri):
    T = regina_tri
    any_progress = False
    progress = True
    while progress:
        progress = False
        edges = sorted(T.edges(), key=lambda e:e.degree())
        for e in edges:
            d = e.degree()
            if d == 1:
                progress = T.twoOneMove(e, 0)
            elif d == 2:
                progress = T.twoZeroMove(e)
            elif d == 3:
                progress = T.threeTwoMove(e)
            else:
                break

            if progress:
                any_progress = True
                break
    return any_progress

def simplify_via_randomization(regina_tri, max_failed_attempts=100):
    T = regina_tri
    any_progress = obvious_simplify(T)
    failures = 0
    while failures < max_failed_attempts:
        edges = [e for e in T.edges() if e.degree() == 4]
        if len(edges) == 0:
            return 
        edge = random.choice(edges)
        T.fourFourMove(edge, random.choice([0, 1]))
        progress = obvious_simplify(T)
        if progress:
            any_progress = True
        else:
            failures += 1

def interesting_two_vertex_triangulation(regina_tri):
    T = regina_tri
    tets = T.tetrahedra()
    for tet in tets:
        T.oneFourMove(tet)

    face = triangles_on_high(T)
    while not face is None:
        T.twoThreeMove(face)
        face = triangles_on_high(T)

    progress = True
    while progress and T.countVertices() > 2:
        progress = False
        edges = sorted(T.edges(), key=lambda e:e.vertex(0).degree() + e.vertex(1).degree())
        for edge in edges:
            if T.collapseEdge(edge):
                progress = True
                break

    obvious_simplify(T)

        
def examine_two_vertex(snappy_manifold):
    M = snappy_manifold
    for iso in util.closed_isosigs(M, 10000, 35):
        for i in range(1):
            R = regina.NTriangulation(iso)
            interesting_two_vertex_triangulation(R)
            new_iso = R.isoSig()
            T = t3m.Mcomplex(new_iso)
            print(len(T))
            for eo in edge_orient.edge_orientations(T):
                if eo.gives_foliation():
                    print(M, new_iso, eo.signs)
                    return 

def basic_examine_two_vertex(iso):
    R = regina.NTriangulation(iso)
    interesting_two_vertex_triangulation(R)
    new_iso = R.isoSig()
    T = t3m.Mcomplex(new_iso)
    print new_iso
    for eo in edge_orient.edge_orientations(T):
        if eo.gives_foliation():
            print(M, new_iso, eo.signs)
            return 

                

if __name__ == '__main__':
    #iso = 'oLLLALAPzQcbedgfjihlkmnlnnxxnxxqaxaxqqnhx'
    #R = regina.NTriangulation(iso)
    #M = snappy.Manifold('m146(6,1)')
    #bad = 'pLALvAzAQzQabcegiklijkmnooobrwbbclhjjnjhxoj_cDbC'
    M = snappy.Manifold('s137(5, 4)')
    #for ex in zhs:
    #    M = snappy.Manifold(ex)
    #    examine_two_vertex(M)

    onevert = 'zLALLzvLLQALzMPQMkcbbeghlpqkpnmtruwsvxxxvyyyqqjnajwahjhcnqrehsqqaqwaet'
    twovert = 'zLLLvwMwMwzPMMLQQkbcgikkllolspttqvrwuyuyxxxytsmjxpiiagpftixavwrakpaxat'

    #Winner for 'm146(6, 1)' sLLLzvQzPPwQQacfglhimkoplnprqrqrnkwkqexxbuavubupipj
    
