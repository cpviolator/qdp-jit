// -*- C++ -*-
// $Id: qdp_scalarsite_sse.h,v 1.12 2003-10-17 20:28:35 edwards Exp $

/*! @file
 * @brief Intel SSE optimizations
 *
 * SSE optimizations of basic operations
 */

#ifndef QDP_SCALARSITE_SSE_H
#define QDP_SCALARSITE_SSE_H

// These SSE asm instructions are only supported under GCC/G++
#if defined(__GNUC__)

QDP_BEGIN_NAMESPACE(QDP);

/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */

// Use this def just to safe some typing later on in the file
#define RComplexFloat  RComplex<float>


typedef struct
{
   unsigned int c1,c2,c3,c4;
} sse_mask __attribute__ ((aligned (16)));

static sse_mask _sse_sgn13 __attribute__ ((unused)) ={0x80000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn24 __attribute__ ((unused)) ={0x00000000, 0x80000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn3  __attribute__ ((unused)) ={0x00000000, 0x00000000, 0x80000000, 0x00000000};
static sse_mask _sse_sgn4  __attribute__ ((unused)) ={0x00000000, 0x00000000, 0x00000000, 0x80000000};


#include "scalarsite_sse/sse_mult_nn.h"
#include "scalarsite_sse/sse_mult_na.h"
#include "scalarsite_sse/sse_mult_an.h"
#include "scalarsite_sse/sse_mat_vec.h"
#include "scalarsite_sse/sse_adj_mat_vec.h"
#include "scalarsite_sse/sse_addvec.h"
#include "scalarsite_sse/sse_mat_hwvec.h"
#include "scalarsite_sse/sse_adj_mat_hwvec.h"

// #define QDP_SCALARSITE_DEBUG

// Optimized version of  
//    PColorMatrix<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * PColorMatrix<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiply>::Type_t
operator*(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	  const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*M" << endl;
#endif

  _inline_sse_mult_su3_nn(l,r,d);

  return d;
}


// Optimized version of  
//    PScalar<PColorMatrix<RComplexFloat,3>> <- PScalar<PColorMatrix<RComplexFloat,3>> * 
//                         PScalar<PColorMatrix<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >,
  PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, OpMultiply>::Type_t
operator*(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	  const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*PSc<M>" << endl;
#endif

  _inline_sse_mult_su3_nn(l.elem(),r.elem(),d.elem());

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- adj(PColorMatrix<RComplexFloat,3>) * PColorMatrix<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	    const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*M" << endl;
#endif

  _inline_sse_mult_su3_an(l,r,d);

  return d;
}


// Optimized version of  
//   PScalar<PColorMatrix<RComplexFloat,3>> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PScalar<PColorMatrix<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	    const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*PSc<M>" << endl;
#endif

  _inline_sse_mult_su3_an(l.elem(),r.elem(),d.elem());

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * adj(PColorMatrix<RComplexFloat,3>)
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiplyAdj>::Type_t
multiplyAdj(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	    const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*adj(M)" << endl;
#endif

  _inline_sse_mult_su3_na(l,r,d);

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * adj(PColorMatrix<RComplexFloat,3>)
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	    const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PS<M>*adj(PSc<M>)" << endl;
#endif

  _inline_sse_mult_su3_na(l.elem(),r.elem(),d.elem());

  return d;
}


#if 0
// Ooops, this macro does not exist!!

// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- adj(PColorMatrix<RComplexFloat,3>) * adj(PColorMatrix<RComplexFloat,3>)
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	       const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiplyAdj>::Type_t  d;

  _inline_sse_mult_su3_aa(l,r,d);

  return d;
}
#endif


// Optimized version of  
//    PColorVector<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * PColorVector<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PVector<RComplexFloat,3,PColorVector>, OpMultiply>::Type_t
operator*(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	  const PVector<RComplexFloat,3,PColorVector>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PVector<RComplexFloat,3,PColorVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*V" << endl;
#endif

  _inline_sse_mult_su3_mat_vec(l,r,d);

  return d;
}


// Optimized version of  
//    PScalar<PColorVector<RComplexFloat,3>> <- PScalar<PColorMatrix<RComplexFloat,3>> * PScalar<PColorVector<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PScalar<PVector<RComplexFloat,3,PColorVector> >, OpMultiply>::Type_t
operator*(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	  const PScalar<PVector<RComplexFloat,3,PColorVector> >& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PScalar<PVector<RComplexFloat,3,PColorVector> >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*PSc<V>" << endl;
#endif

  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(),d.elem());

  return d;
}


// Optimized version of  
//    PColorVector<RComplexFloat,3> <- adj(PColorMatrix<RComplexFloat,3>) * PColorVector<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PVector<RComplexFloat,3,PColorVector>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	    const PVector<RComplexFloat,3,PColorVector>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PVector<RComplexFloat,3,PColorVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*V" << endl;
#endif

  _inline_sse_mult_adj_su3_mat_vec(l,r,d);

  return d;
}


// Optimized version of   StaggeredFermion <- ColorMatrix*StaggeredFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,1> <- PScalar<PColorMatrix<RComplexFloat,3>> * PSpinVector<PColorVector<RComplexFloat,3>,1>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PVector<PVector<RComplexFloat,3,PColorVector>,1,PSpinVector>, OpMultiply>::Type_t
operator*(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	  const PVector<PVector<RComplexFloat,3,PColorVector>,1,PSpinVector>& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PVector<PVector<RComplexFloat,3,PColorVector>,1,PSpinVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*S" << endl;
#endif

  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));

  return d;
}


// Optimized version of  
//    PScalar<PColorVector<RComplexFloat,3>> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PScalar<PColorVector<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PScalar<PVector<RComplexFloat,3,PColorVector> >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	    const PScalar<PVector<RComplexFloat,3,PColorVector> >& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PScalar<PVector<RComplexFloat,3,PColorVector> >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*PSc<V>" << endl;
#endif

  _inline_sse_mult_adj_su3_mat_vec(l.elem(),r.elem(),d.elem());

  return d;
}


// Optimized version of   StaggeredFermion <- adj(ColorMatrix)*StaggeredFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,1> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PSpinVector<PColorVector<RComplexFloat,3>,1>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PVector<PVector<RComplexFloat,3,PColorVector>,1,PSpinVector>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	    const PVector<PVector<RComplexFloat,3,PColorVector>,1,PSpinVector>& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PVector<PVector<RComplexFloat,3,PColorVector>,1,PSpinVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(S)*S" << endl;
#endif

  _inline_sse_mult_adj_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,2> <- PScalar<PColorMatrix<RComplexFloat,3>> * 
//                     PSpinVector<ColorVector<RComplexFloat,3>,2>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PVector<PVector<RComplexFloat,3,PColorVector>,2,PSpinVector>, OpMultiply>::Type_t
operator*(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
          const PVector<PVector<RComplexFloat,3,PColorVector>,2,PSpinVector>& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PVector<PVector<RComplexFloat,3,PColorVector>,2,PSpinVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*H" << endl;
#endif

#if 0
  _inline_sse_mult_su3_mat_hwvec(l,r,d);
#else
  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));
#endif

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,2> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * 
//                     PSpinVector<ColorVector<RComplexFloat,3>,2>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PVector<PVector<RComplexFloat,3,PColorVector>,2,PSpinVector>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
            const PVector<PVector<RComplexFloat,3,PColorVector>,2,PSpinVector>& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PVector<PVector<RComplexFloat,3,PColorVector>,2,PSpinVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*H" << endl;
#endif

  _inline_sse_mult_adj_su3_mat_hwvec(l,r,d);

  return d;
}


// Optimized version of  
//    PColorVector<RComplexFloat,3> <- PColorVector<RComplexFloat,3> + PColorVector<RComplexFloat,3>
template<>
inline BinaryReturn<PVector<RComplexFloat,3,PColorVector>, 
  PVector<RComplexFloat,3,PColorVector>, OpAdd>::Type_t
operator+(const PVector<RComplexFloat,3,PColorVector>& l, 
	  const PVector<RComplexFloat,3,PColorVector>& r)
{
  BinaryReturn<PVector<RComplexFloat,3,PColorVector>, 
    PVector<RComplexFloat,3,PColorVector>, OpAdd>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "V+V" << endl;
#endif

  _inline_sse_add_su3_vector(l,r,d);

  return d;
}


// Optimized version of   DiracFermion <- ColorMatrix*DiracFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,4> <- PScalar<PColorMatrix<RComplexFloat,3>> 
//                           * PSpinVector<PColorVector<RComplexFloat,3>,4>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PVector<PVector<RComplexFloat,3,PColorVector>,4,PSpinVector>, OpMultiply>::Type_t
operator*(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	  const PVector<PVector<RComplexFloat,3,PColorVector>,4,PSpinVector>& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PVector<PVector<RComplexFloat,3,PColorVector>,4,PSpinVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*D" << endl;
#endif

  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));
  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(2),d.elem(2));
  _inline_sse_mult_su3_mat_vec(l.elem(),r.elem(3),d.elem(3));

  return d;
}


// Optimized version of   DiracFermion <- adj(ColorMatrix)*DiracFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,4> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>)
//                           * PSpinVector<PColorVector<RComplexFloat,3>,4>
template<>
inline BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
  PVector<PVector<RComplexFloat,3,PColorVector>,4,PSpinVector>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >& l, 
	    const PVector<PVector<RComplexFloat,3,PColorVector>,4,PSpinVector>& r)
{
  BinaryReturn<PScalar<PMatrix<RComplexFloat,3,PColorMatrix> >, 
    PVector<PVector<RComplexFloat,3,PColorVector>,4,PSpinVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*D" << endl;
#endif

  _inline_sse_mult_adj_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_sse_mult_adj_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));
  _inline_sse_mult_adj_su3_mat_vec(l.elem(),r.elem(2),d.elem(2));
  _inline_sse_mult_adj_su3_mat_vec(l.elem(),r.elem(3),d.elem(3));

  return d;
}



// Optimized version of  
//    PScalar<PColorVector<RComplexFloat,3>> <- PScalar<PColorVector<RComplexFloat,3>> 
//                                            + PScalar<PColorVector<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PVector<RComplexFloat,3,PColorVector> >, 
  PScalar<PVector<RComplexFloat,3,PColorVector> >, OpAdd>::Type_t
operator+(const PScalar<PVector<RComplexFloat,3,PColorVector> >& l, 
	  const PScalar<PVector<RComplexFloat,3,PColorVector> >& r)
{
  BinaryReturn<PScalar<PVector<RComplexFloat,3,PColorVector> >, 
    PScalar<PVector<RComplexFloat,3,PColorVector> >, OpAdd>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<V>+PSc<V>" << endl;
#endif

  _inline_sse_add_su3_vector(l.elem(),r.elem(),d.elem());

  return d;
}


#if 1
// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * LatticeColorMatrix
// NOTE: let this be a subroutine to save space
template<>
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const OrderedSubset& s);

#endif

#if 1
// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * LatticeColorMatrix
// NOTE: let this be a subroutine to save space
template<>
inline 
void evaluate(OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PSpinVector<PColorVector<RComplexFloat, 3>, 2>, 
	      OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > > > >,
	      OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > >& rhs,
	      const OrderedSubset& s)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "specialized QDP_H_M_times_H" << endl;
#endif

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  typedef OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > H;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const H& r = static_cast<const H&>(rhs.expression().right());

  for(int i=s.start(); i <= s.end(); ++i) 
  {
#if 0
    // This form appears significantly slower than below
    _inline_sse_mult_su3_mat_hwvec(l.elem(i),
				   r.elem(i),
				   d.elem(i));
#else
    _inline_sse_mult_su3_mat_vec(l.elem(i).elem(),
				 r.elem(i).elem(0),
				 d.elem(i).elem(0));
    _inline_sse_mult_su3_mat_vec(l.elem(i).elem(),
				 r.elem(i).elem(1),
				 d.elem(i).elem(1));
#endif
  }
}

#endif

/*! @} */   // end of group optimizations

#if defined(QDP_SCALARSITE_DEBUG)
#undef QDP_SCALARSITE_DEBUG
#endif

QDP_END_NAMESPACE();

#endif  // defined(__GNUC__)

#endif
