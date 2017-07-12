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
    for iso in util.closed_isosigs(M)[:10]:
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
    M = snappy.Manifold('m231(-6, 1)')
    examples = ['m146(6, 1)', 'm231(-6, 1)', 'm232(-7, 1)', 'm233(-6, 1)', 'm241(1, 5)', 'm243(-4, 3)', 'm249(6, 1)', 'm253(1, 4)', 'm257(1, 5)', 'm280(1, 5)', 'm290(4, 1)', 'm290(5, 1)', 'm318(-1, 5)', 'm336(-1, 5)', 'm337(-5, 1)', 'm337(-6, 1)', 'm349(-5, 1)', 'm349(-6, 1)', 'm351(-4, 5)', 'm360(-1, 5)']
    
    known = ['m003(-1, 3)', 'm003(-1, 4)', 'm003(-1, 5)', 'm003(-2, 5)', 'm004(-1, 2)', 'm004(-3, 2)', 'm004(-5, 2)', 'm004(-6, 1)', 'm004(-7, 1)', 'm004(-7, 2)', 'm004(-8, 1)', 'm004(-9, 1)', 'm006(-5, 1)', 'm006(-5, 2)', 'm006(-6, 1)', 'm006(-7, 1)', 'm007(-5, 1)', 'm007(-5, 2)', 'm007(-6, 1)', 'm007(-7, 1)']

    zhs = ['m232(-7, 1)', 'm395(-5, 1)', 's517(-1, 4)', 's547(3, 4)', 's636(-4, 3)', 's642(-4, 5)', 's668(-6, 1)', 's805(4, 5)', 'v1181(7, 1)', 'v1828(6, 1)', 'v1925(4, 3)', 'v1931(5, 1)', 'v2397(5, 4)', 'v2427(2, 5)', 'v2841(5, 1)', 'v2886(-5, 1)', 'v2993(-6, 1)', 'v3006(-4, 3)', 'v3177(-5, 1)', 'v3307(-1, 5)', 'v3484(-5, 1)', 't03313(-7, 1)', 't05001(5, 2)', 't05144(5, 1)', 't05259(5, 1)', 't06908(5, 3)', 't07039(-5, 1)', 't07218(-6, 1)', 't07311(-3, 4)', 't07728(-7, 1)', 't07873(-6, 1)', 't08142(-5, 1)', 't08210(-5, 1)', 't08260(5, 3)', 't08653(4, 3)', 't08833(-4, 5)', 't08905(-4, 1)', 't08968(-6, 1)', 't09044(8, 1)', 't09097(1, 5)', 't09116(4, 5)', 't09349(-4, 1)', 't09627(-5, 1)', 't09779(-4, 1)', 't09854(3, 5)', 't10481(8, 1)', 't10482(1, 4)', 't10833(-5, 1)', 't10966(-1, 5)', 't10990(-5, 1)', 't11075(-5, 2)', 't11545(3, 1)', 't12134(-4, 3)', 't12196(6, 1)', 't12278(5, 1)', 't12575(-7, 1)', 'o9_11006(1, 5)', 'o9_16582(-5, 2)', 'o9_21144(5, 2)', 'o9_21970(7, 1)', 'o9_22265(-4, 1)', 'o9_22322(1, 4)', 'o9_23232(7, 1)', 'o9_23541(5, 4)', 'o9_24096(2, 5)', 'o9_24232(-2, 3)', 'o9_24291(-4, 3)', 'o9_25033(-4, 5)', 'o9_25247(8, 1)', 'o9_25545(3, 4)', 'o9_25791(6, 1)', 'o9_25900(4, 3)', 'o9_26816(6, 1)', 'o9_26867(-5, 1)', 'o9_26889(1, 5)', 'o9_26948(5, 4)', 'o9_28029(4, 1)', 'o9_28479(1, 4)', 'o9_28546(-2, 3)', 'o9_28726(5, 4)', 'o9_28965(4, 1)', 'o9_29301(-5, 1)', 'o9_29337(4, 1)', 'o9_29573(7, 1)', 'o9_29756(-4, 3)', 'o9_29825(6, 1)', 'o9_30056(5, 1)', 'o9_30119(2, 3)', 'o9_30463(-4, 3)', 'o9_30560(1, 5)', 'o9_30775(-4, 1)', 'o9_32093(4, 3)', 'o9_32482(-5, 4)', 'o9_32569(1, 3)', 'o9_32864(4, 1)', 'o9_32893(-8, 1)', 'o9_33014(5, 2)', 'o9_33595(-2, 3)', 'o9_33712(5, 1)', 'o9_34835(4, 5)', 'o9_34852(1, 5)', 'o9_35339(5, 1)', 'o9_35390(3, 2)', 'o9_35479(-3, 4)', 'o9_35883(5, 1)', 'o9_35950(4, 1)', 'o9_36046(4, 1)', 'o9_36352(-5, 1)', 'o9_36379(-5, 1)', 'o9_36984(5, 1)', 'o9_37192(1, 5)', 'o9_37215(3, 2)', 'o9_37360(6, 1)', 'o9_37457(6, 1)', 'o9_38215(2, 3)', 'o9_38268(-1, 5)', 'o9_38532(6, 1)', 'o9_38693(-2, 3)', 'o9_38799(-4, 3)', 'o9_39482(3, 5)', 'o9_39511(6, 1)', 'o9_39656(3, 4)', 'o9_39985(1, 3)', 'o9_40491(-1, 5)', 'o9_40829(1, 3)', 'o9_40940(-5, 1)', 'o9_41042(5, 2)', 'o9_41185(-4, 1)', 'o9_41436(-5, 1)', 'o9_41461(5, 1)', 'o9_41687(-1, 3)', 'o9_41746(-4, 1)', 'o9_41850(5, 1)', 'o9_41896(-5, 1)', 'o9_41967(-4, 1)', 'o9_42157(-2, 3)', 'o9_42570(-2, 3)', 'o9_42571(-3, 5)', 'o9_42582(5, 2)', 'o9_42585(-1, 3)', 'o9_43091(4, 3)', 'o9_43177(-6, 1)', 'o9_43568(1, 3)', 'o9_43648(-4, 3)', 'o9_43793(-1, 3)']

    #for ex in zhs:
    #    M = snappy.Manifold(ex)
    #    examine_two_vertex(M)

    onevert = 'zLALLzvLLQALzMPQMkcbbeghlpqkpnmtruwsvxxxvyyyqqjnajwahjhcnqrehsqqaqwaet'
    twovert = 'zLLLvwMwMwzPMMLQQkbcgikkllolspttqvrwuyuyxxxytsmjxpiiagpftixavwrakpaxat'

    #Winner for 'm146(6, 1)' sLLLzvQzPPwQQacfglhimkoplnprqrqrnkwkqexxbuavubupipj
    
