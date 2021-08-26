// -*- C++ -*-

/*! \file
 * \brief Primitive classes
 *
 * Primitives are the various types on the fibers at the lattice sites
 */


#ifndef QDP_PRIMITIVE_H
#define QDP_PRIMITIVE_H


/*! \defgroup fiber Fiber only types and operations
 * \ingroup fiberbundle
 *
 * Primitives are the various types on the fibers at the lattice sites.
 *
 * The primitive indices, including Reality (also known as complex or real),
 * is represented as a tensor product over various vector spaces. Different
 * kinds of object can transform in those vector spaces, like Scalar, Vector, and
 * Matrix.
 */

#include "qdp_primscalarjit.h"
#include "qdp_primscalarreg.h"
#include "qdp_primscalar.h"

#include "qdp_primmatrixjit.h"
#include "qdp_primmatrixreg.h"
#include "qdp_primmatrix.h"

#include "qdp_primvectorjit.h"
#include "qdp_primvectorreg.h"
#include "qdp_primvector.h"

#include "qdp_primseedjit.h"
#include "qdp_primseedreg.h"
#include "qdp_primseed.h"

#include "qdp_primcolormatjit.h"
#include "qdp_primcolormatreg.h"
#include "qdp_primcolormat.h"

#include "qdp_primcolorvecjit.h"
#include "qdp_primcolorvecreg.h"
#include "qdp_primcolorvec.h"

#include "qdp_primspinmatjit.h"
#include "qdp_primspinmatreg.h"
#include "qdp_primspinmat.h"

#include "qdp_primspinvecjit.h"
#include "qdp_primspinvecreg.h"
#include "qdp_primspinvec.h"

#endif
