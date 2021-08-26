// -*- C++ -*-

/*! \file
 * \brief Primitive Spin Matrix
 */

#ifndef QDP_PRIMSPINMATREG_H
#define QDP_PRIMSPINMATREG_H

namespace QDP {


//-------------------------------------------------------------------------------------
/*! \addtogroup primspinmatrix Spin matrix primitive
 * \ingroup primmatrix
 *
 * Primitive type that transforms like a Spin Matrix
 *
 * @{
 */


//! Primitive Spin Matrix class
template <class T, int N> class PSpinMatrixREG : public PMatrixREG<T, N, PSpinMatrixREG>
{
  //  PSpinMatrixREG(const PSpinMatrixREG& a);
public:
  PSpinMatrixREG(){}

  void setup( const typename JITType< PSpinMatrixREG >::Type_t& j ) {
    for (int i = 0 ; i < N ; i++ ) 
      for (int q = 0 ; q < N ; q++ ) 
	this->elem(i,q).setup( j.elem(i,q) );
  }

  void setup_value( const typename JITType< PSpinMatrixREG >::Type_t& j ) {
    for (int i = 0 ; i < N ; i++ ) 
      for (int q = 0 ; q < N ; q++ ) 
	this->elem(i,q).setup_value( j.elem(i,q) );
  }


  template<class T1>
  PSpinMatrixREG(const PSpinMatrixREG<T1,N>& a)
  {
    this->assign(a);
  }

  //! PSpinMatrixREG = PScalarREG
  /*! Fill with primitive scalar */
  template<class T1>
  inline
  PSpinMatrixREG& operator=(const PScalarREG<T1>& rhs)
    {
      this->assign(rhs);
      return *this;
    }

  //! PSpinMatrixREG = PSpinMatrixREG
  /*! Set equal to another PSpinMatrixREG */
  template<class T1>
  inline
  PSpinMatrixREG& operator=(const PSpinMatrixREG<T1,N>& rhs) 
    {
      this->assign(rhs);
      return *this;
    }


  PSpinMatrixREG& operator=(const PSpinMatrixREG& rhs) 
    {
      this->assign(rhs);
      return *this;
    }

};

/*! @} */   // end of group primspinmatrix




//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

template<class T1, int N>
struct JITType<PSpinMatrixREG<T1,N> > 
{
  typedef PSpinMatrixJIT<typename JITType<T1>::Type_t,N>  Type_t;
};



// Underlying word type
template<class T1, int N>
struct WordType<PSpinMatrixREG<T1,N> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

template<class T1, int N>
struct SinglePrecType<PSpinMatrixREG<T1, N> >
{
  typedef PSpinMatrixREG< typename SinglePrecType<T1>::Type_t , N > Type_t;
};

template<class T1, int N>
struct DoublePrecType<PSpinMatrixREG<T1, N> >
{
  typedef PSpinMatrixREG< typename DoublePrecType<T1>::Type_t , N > Type_t;
};

// Internally used scalars
template<class T, int N>
struct InternalScalar<PSpinMatrixREG<T,N> > {
  typedef PScalarREG<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive into a scalar leaving grid alone
template<class T, int N>
struct PrimitiveScalar<PSpinMatrixREG<T,N> > {
  typedef PScalarREG<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T, int N>
struct LatticeScalar<PSpinMatrixREG<T,N> > {
  typedef PSpinMatrixREG<typename LatticeScalar<T>::Type_t, N>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PSpinMatrixREG) -> PSpinMatrixREG
template<class T1, int N, class Op>
struct UnaryReturn<PSpinMatrixREG<T1,N>, Op> {
  typedef PSpinMatrixREG<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};

// Default binary(PScalarREG,PSpinMatrixREG) -> PSpinMatrixREG
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, Op> {
  typedef PSpinMatrixREG<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinMatrixREG,PSpinMatrixREG) -> PSpinMatrixREG
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, Op> {
  typedef PSpinMatrixREG<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinMatrixREG,PScalarREG) -> PSpinMatrixREG
template<class T1, int N, class T2, class Op>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, Op> {
  typedef PSpinMatrixREG<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<PSpinMatrixREG<T2,N>, OpCast<T1> > {
  typedef PScalarREG<typename UnaryReturn<T, OpCast>::Type_t, N>  Type_t;
//  typedef T1 Type_t;
};
#endif


// Assignment is different
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, OpAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, OpAddAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, OpSubtractAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, OpMultiplyAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, OpAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, OpAddAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, OpSubtractAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, OpMultiplyAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, OpDivideAssign > {
  typedef PSpinMatrixREG<T1,N> &Type_t;
};
 


// SpinMatrix
template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnTrace > {
  typedef PScalarREG<typename UnaryReturn<T, FnTrace>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnRealTrace > {
  typedef PScalarREG<typename UnaryReturn<T, FnRealTrace>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnImagTrace > {
  typedef PScalarREG<typename UnaryReturn<T, FnImagTrace>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnNorm2 > {
  typedef PScalarREG<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnLocalNorm2 > {
  typedef PScalarREG<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnTraceMultiply> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnTraceMultiply> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnTraceMultiply> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnInnerProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnInnerProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnInnerProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnLocalInnerProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};



template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnLocalColorInnerProduct>::Type_t
localColorInnerProduct(const PSpinMatrixREG<T1,N>& l, const PSpinMatrixREG<T2,N>& r)
{
  typename BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnLocalColorInnerProduct>::Type_t d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      {
	d.elem(i,j) = localColorInnerProduct(l.elem(0,i), r.elem(0,j));
	for(int k=1; k < N; ++k)
	  d.elem(i,j) += localColorInnerProduct(l.elem(k,i), r.elem(k,j));
      }
    
  return d;
}


  
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnLocalInnerProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnLocalInnerProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnInnerProductReal> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnInnerProductReal> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnInnerProductReal> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnLocalInnerProductReal> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnLocalInnerProductReal> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnLocalInnerProductReal> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};






//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------
/*! \addtogroup primspinmatrix */
/*! @{ */

// SpinMatrix class primitive operations

// trace = traceSpin(source1)
/*! This only acts on spin indices and is diagonal in all other indices */
template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnTraceSpin > {
  typedef PScalarREG<typename UnaryReturn<T, FnTraceSpin>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PSpinMatrixREG<T,N>, FnTraceSpin>::Type_t
traceSpin(const PSpinMatrixREG<T,N>& s1)
{
  typename UnaryReturn<PSpinMatrixREG<T,N>, FnTraceSpin>::Type_t  d;
  
  // Since the spin index is eaten, do not need to pass on function by
  // calling trace(...) again
  d.elem() = s1.elem(0,0);
  for(int i=1; i < N; ++i)
    d.elem() += s1.elem(i,i);

  return d;
}

//! traceSpinMultiply(source1,source2)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnTraceSpinMultiply> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PSpinMatrixREG<T1,N>& l, const PSpinMatrixREG<T2,N>& r)
{
  typename BinaryReturn<PSpinMatrixREG<T1,N>, PSpinMatrixREG<T2,N>, FnTraceSpinMultiply>::Type_t  d;

  // The traceSpin is eaten here
  d.elem() = l.elem(0,0) * r.elem(0,0);
  for(int k=1; k < N; ++k)
    d.elem() += l.elem(0,k) * r.elem(k,0);

  for(int j=1; j < N; ++j)
    for(int k=0; k < N; ++k)
      d.elem() += l.elem(j,k) * r.elem(k,j);

  return d;
}

//! PScalarREG = traceSpinMultiply(PSpinMatrixREG,PScalarREG)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnTraceSpinMultiply> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PSpinMatrixREG<T1,N>& l, const PScalarREG<T2>& r)
{
  typename BinaryReturn<PSpinMatrixREG<T1,N>, PScalarREG<T2>, FnTraceSpinMultiply>::Type_t  d;

  // The traceSpin is eaten here
  d.elem() = l.elem(0,0) * r.elem();
  for(int k=1; k < N; ++k)
    d.elem() += l.elem(k,k) * r.elem();

  return d;
}

// PScalarREG = traceSpinMultiply(PScalarREG,PSpinMatrixREG)
template<class T1, class T2, int N>
struct BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnTraceSpinMultiply> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PScalarREG<T1>& l, const PSpinMatrixREG<T2,N>& r)
{
  typename BinaryReturn<PScalarREG<T1>, PSpinMatrixREG<T2,N>, FnTraceSpinMultiply>::Type_t  d;

  // The traceSpin is eaten here
  d.elem() = l.elem() * r.elem(0,0);
  for(int k=1; k < N; ++k)
    d.elem() += l.elem() * r.elem(k,k);

  return d;
}



/*! Specialise the return type */
template <class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnTransposeSpin > {
  typedef PSpinMatrixREG<typename UnaryReturn<T, FnTransposeSpin>::Type_t, N> Type_t;
};

//! PSpinMatrixREG = transposeSpin(PSpinMatrixREG) 
/*! t = transposeSpin(source1) - SpinMatrix specialization -- where the work is actually done */
template<class T, int N>
inline typename UnaryReturn<PSpinMatrixREG<T,N>, FnTransposeSpin >::Type_t
transposeSpin(const PSpinMatrixREG<T,N>& s1)
{
  typename UnaryReturn<PSpinMatrixREG<T,N>, FnTransposeSpin>::Type_t d;
 
  for(int i=0; i < N; i++) { 
    for(int j=0; j < N; j++) { 
      // Transpose, so flip indices
      d.elem(i,j) = s1.elem(j,i);
    }
  }
  return d;
}


//-----------------------------------------------
// OuterProduct must be handled specially for each color and spin
// The problem is the traits class - I have no way to say to PVectorREG's
//  transform into a PMatrixREG but downcast the trait to a PColorMatrix or
//  PSpinMatrixREG

//! PSpinMatrixREG = outerProduct(PSpinVectorREG, PSpinVectorREG)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVectorREG<T1,N>, PSpinVectorREG<T2,N>, FnOuterProduct> {
  typedef PSpinMatrixREG<typename BinaryReturn<T1, T2, FnOuterProduct>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVectorREG<T1,N>, PSpinVectorREG<T2,N>, FnOuterProduct>::Type_t
outerProduct(const PSpinVectorREG<T1,N>& l, const PSpinVectorREG<T2,N>& r)
{
  typename BinaryReturn<PSpinVectorREG<T1,N>, PSpinVectorREG<T2,N>, FnOuterProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = outerProduct(l.elem(i),r.elem(j));

  return d;
}


//-----------------------------------------------
// Optimization of traceSpin(outerProduct(PSpinVectorREG, PSpinVectorREG))

//! PScalarREG = traceSpinOuterProduct(PSpinVectorREG, PSpinVectorREG)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVectorREG<T1,N>, PSpinVectorREG<T2,N>, FnTraceSpinOuterProduct> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnOuterProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVectorREG<T1,N>, PSpinVectorREG<T2,N>, FnTraceSpinOuterProduct>::Type_t
traceSpinOuterProduct(const PSpinVectorREG<T1,N>& l, const PSpinVectorREG<T2,N>& r)
{
  typename BinaryReturn<PSpinVectorREG<T1,N>, PSpinVectorREG<T2,N>, FnTraceSpinOuterProduct>::Type_t  d;

  d.elem() = outerProduct(l.elem(0),r.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += outerProduct(l.elem(i),r.elem(i));

  return d;
}


//-----------------------------------------------
// Peeking and poking
//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T, int N>
struct UnaryReturn<PSpinMatrixREG<T,N>, FnPeekSpinMatrixREG > {
  typedef PScalarREG<typename UnaryReturn<T, FnPeekSpinMatrixREG >::Type_t>  Type_t;
};


template<class T, int N>
inline typename UnaryReturn<PSpinMatrixREG<T,N>, FnPeekSpinMatrixREG >::Type_t
peekSpin(const PSpinMatrixREG<T,N>& l, llvm::Value* row, llvm::Value* col)
{
  typename UnaryReturn<PSpinMatrixREG<T,N>, FnPeekSpinMatrixREG >::Type_t  d;

  typedef typename JITType< PSpinMatrixREG<T,N> >::Type_t TTjit;

  llvm::Value* ptr_local = llvm_alloca( llvm_get_type<typename WordType<T>::Type_t>() , TTjit::Size_t );

  TTjit dj;
  dj.setup( ptr_local, JitDeviceLayout::Scalar );
  dj=l;

  d.elem() = dj.getRegElem(row,col);
  return d;
}

//! Insert spin matrix components
template<class T1, class T2, int N>
inline PSpinMatrixREG<T1,N>&
pokeSpin(PSpinMatrixREG<T1,N>& l, const PScalarREG<T2>& r, int row, int col)
{
  // Note, do not need to propagate down since the function is eaten at this level
  l.getRegElem(row,col) = r.elem();
  return l;
}



//-----------------------------------------------------------------------------
//! PSpinVectorREG<T,4> = P_+ * PSpinVectorREG<T,4>
template<class T>
inline typename UnaryReturn<PSpinMatrixREG<T,4>, FnChiralProjectPlus>::Type_t
chiralProjectPlus(const PSpinMatrixREG<T,4>& s1)
{
  typename UnaryReturn<PSpinMatrixREG<T,4>, FnChiralProjectPlus>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = s1.elem(0,i);
    d.elem(1,i) = s1.elem(1,i);
    zero_rep(d.elem(2,i));
    zero_rep(d.elem(3,i));
  }

  return d;
}

//! PSpinVectorREG<T,4> = P_- * PSpinVectorREG<T,4>
template<class T>
inline typename UnaryReturn<PSpinMatrixREG<T,4>, FnChiralProjectMinus>::Type_t
chiralProjectMinus(const PSpinMatrixREG<T,4>& s1)
{
  typename UnaryReturn<PSpinMatrixREG<T,4>, FnChiralProjectMinus>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    zero_rep(d.elem(0,i));
    zero_rep(d.elem(1,i));
    d.elem(2,i) = s1.elem(2,i);
    d.elem(3,i) = s1.elem(3,i);
  }

  return d;
}

//------------------------------------------
// PScalarREG = traceSpinQuarkContract13(PSpinMatrixREG,PSpinMatrixREG)
template<class T1, class T2>
struct BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnTraceSpinQuarkContract13> {
  typedef PScalarREG<typename BinaryReturn<T1, T2, FnTraceSpinQuarkContract13>::Type_t>  Type_t;
};

//! PScalarREG = traceSpinQuarkContract13(PSpinMatrixREG,PSpinMatrixREG)
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnTraceSpinQuarkContract13>::Type_t
traceSpinQuarkContract13(const PSpinMatrixREG<T1,4>& l, const PSpinMatrixREG<T2,4>& r)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnTraceSpinQuarkContract13>::Type_t  d;

  d.elem() = quarkContractXX(l.elem(0,0), r.elem(0,0));
  for(int k=1; k < 4; ++k)
    d.elem() += quarkContractXX(l.elem(k,0), r.elem(k,0));

  for(int j=1; j < 4; ++j)
    for(int k=0; k < 4; ++k)
      d.elem() += quarkContractXX(l.elem(k,j), r.elem(k,j));

  return d;
}


// quark propagator contraction
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract13>::Type_t
quarkContract13(const PSpinMatrixREG<T1,4>& s1, const PSpinMatrixREG<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract13>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(0,i), s2.elem(0,j));
      for(int k=1; k < 4; ++k)
       	d.elem(i,j) += quarkContractXX(s1.elem(k,i), s2.elem(k,j));
    }
  
  return d;
}
  

template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract14>::Type_t
quarkContract14(const PSpinMatrixREG<T1,4>& s1, const PSpinMatrixREG<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract14>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(0,i), s2.elem(j,0));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(k,i), s2.elem(j,k));
    }

  return d;
}

  
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract23>::Type_t
quarkContract23(const PSpinMatrixREG<T1,4>& s1, const PSpinMatrixREG<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract23>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(i,0), s2.elem(0,j));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(i,k), s2.elem(k,j));
    }
  
  return d;
}

  
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract24>::Type_t
quarkContract24(const PSpinMatrixREG<T1,4>& s1, const PSpinMatrixREG<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract24>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(i,0), s2.elem(j,0));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(i,k), s2.elem(j,k));
    }
  
  return d;
}

  
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract12>::Type_t
quarkContract12(const PSpinMatrixREG<T1,4>& s1, const PSpinMatrixREG<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract12>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(0,0), s2.elem(i,j));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(k,k), s2.elem(i,j));
    }

  return d;
}

  
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract34>::Type_t
quarkContract34(const PSpinMatrixREG<T1,4>& s1, const PSpinMatrixREG<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrixREG<T1,4>, PSpinMatrixREG<T2,4>, FnQuarkContract34>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(i,j), s2.elem(0,0));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(i,j), s2.elem(k,k));
    }

  return d;
}

/*! @} */   // end of group primspinmatrix

} // namespace QDP

#endif
