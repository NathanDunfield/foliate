import pyregina
import snappy
M = snappy.Manifold('m004')
print('Regina version: %s' % pyregina.version())
R = pyregina.Triangulation(M)
print(R.num_tetrahedra())
print(R.isosig())
print(R.retriangulate(5, True))





