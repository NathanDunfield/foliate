Foliar orientations
===================

To install the "foliar" module into SageMath+SnapPy do the following
in this directory::

  sage -pip install .

You can test to make sure all is well by doing::

  sage -python -m foliar.test

which should conclude with a message such as::

  All doctests:
    0 failures out of 73 tests.

Example of searching for a foliar orientation on various
triangulations of a closed manifold::

  sage: import snappy, foliar
  sage: M = snappy.Manifold('m004(1,2)')
  sage: eo = foliar.first_foliation(M, 5, 25)
  sage: eo.gives_foliation()
  True
  sage: eo.euler_class_vanishes()
  True

Example of searching for a persistently foliar orientation on a
1-cusped manifold::

  sage: m004 = snappy.Manifold('m004')
  sage: foliar.degeneracy_slopes(m004)
  [(1, 0)]
  sage: m003 = snappy.Manifold('m003')
  sage: foliar.degeneracy_slopes_with_search(m003)
  ([], [])

Note here that m003 is Floer simple and hence nothing is found.

For many more useable examples, see the docstrings the files
"src/main.py" and "src/edge_orient.py".

