// $Id: qdp_scalarsite_sse.cc,v 1.5 2003-10-22 02:44:46 edwards Exp $

/*! @file
 * @brief Intel SSE optimizations
 * 
 * SSE optimizations of basic operations
 */


#include "qdp.h"

// These SSE asm instructions are only supported under GCC/G++
#if defined(__GNUC__)

QDP_BEGIN_NAMESPACE(QDP);

// Specialization to optimize the case   
//    LatticeColorMatrix[OrderedSubset] = LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > >, 
	      OLattice<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      Reference<QDPType<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > >, 
	      OLattice<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > > > > > >,
	      OLattice<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > > > >& rhs,
	      const OrderedSubset& s)
{
//  cout << "call single site QDP_M_eq_M_times_M" << endl;

  typedef OLattice<PScalar<PScalar<PColorMatrix<RComplexFloat, 3> > > >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  for(int i=s.start(); i <= s.end(); ++i) 
  {
    _inline_sse_mult_su3_nn(l.elem(i).elem().elem(),r.elem(i).elem().elem(),d.elem(i).elem().elem());
  }
}

QDP_END_NAMESPACE();

#endif  // defined(__GNUC__)
