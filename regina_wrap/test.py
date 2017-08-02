import pyregina
import snappy
M = snappy.Manifold('m004')
print pyregina.version()
M = pyregina.Triangulation(M._to_string())
print M.num_tetrahedra()
print M.isosig()
M.retriangulate(6)
#print n
#print M.fundamental_group()




