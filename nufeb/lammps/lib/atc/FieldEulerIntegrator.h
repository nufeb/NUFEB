#ifndef FIELD_EULER_INTEGRATOR_H
#define FIELD_EULER_INTEGRATOR_H

// ATC includes
#include "Array2D.h"
#include "MatrixLibrary.h"
#include "PhysicsModel.h"
#include "TimeIntegrator.h"
#include "ImplicitSolveOperator.h"

// other includes
#include <vector>
#include <map>

namespace ATC {

// Forward class declarations
class ATC_Coupling;
class FE_Engine;

/**
 *  @class FieldEulerIntegrator
 *  @brief method for integrating fast fields
 */

class FieldEulerIntegrator {

 public:

  /** Constructor */
  FieldEulerIntegrator(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    /*const*/ FE_Engine * feEngine,
    /*const*/ ATC_Coupling * atc,
    const Array2D< bool > & rhsMask  // copy
  );


  /** Destructor */
  virtual ~FieldEulerIntegrator() {};

  /** initialize */
  virtual void initialize(const double /* dt */, const double /* time */,
                          FIELDS & /* fields */) {};

  /** update */
  virtual void update(const double dt, const double time,
    FIELDS & fields,  FIELDS & rhs) = 0;

 protected:

  /** Pointer to ATC_Tranfer */
  ATC_Coupling * atc_;

  /** Pointer to FE_Engine */
  /*const*/ FE_Engine * feEngine_;

  /** Pointer to PhysicsModel */
  const PhysicsModel * physicsModel_;

  /** field name */
  FieldName fieldName_; 

  /** rhs mask */
  Array2D <bool> rhsMask_;

  /** number of nodes */
  int nNodes_;
};

/**
 *  @class FieldExplicitEulerIntegrator
 *  @brief explicit Euler method for integrating fast electron fields
 */
class FieldExplicitEulerIntegrator : public FieldEulerIntegrator {

 public:

  /** Constructor */
  FieldExplicitEulerIntegrator(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    /*const*/ FE_Engine * feEngine,
    /*const*/ ATC_Coupling * atc,
    const Array2D< bool > & rhsMask  // copy
  );

  /** Destructor */
  virtual ~FieldExplicitEulerIntegrator() {};

  /** update */
  void update(const double dt, const double time,
    FIELDS & fields,  FIELDS & rhs);

};

/**
 *  @class FieldImplicitEulerIntegrator
 *  @brief explicit Euler method for integrating fast electron fields
 */
class FieldImplicitEulerIntegrator : public FieldEulerIntegrator {

 public:

  /** Constructor */
  FieldImplicitEulerIntegrator(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    /*const*/ FE_Engine * feEngine,
    /*const*/ ATC_Coupling * atc,
    const Array2D< bool > & rhsMask, // copy
    const double alpha = 0.5 // default to trap/midpt
  );

  /** Destructor */
  virtual ~FieldImplicitEulerIntegrator() {};

  /** update */
  void update(const double dt, const double time,
    FIELDS & fields,  FIELDS & rhs);

 protected:
  /** euler update factor */
  double alpha_;

  /** perturbation */
  double dT_;

  /** max number of restarts = size of basis */
  int maxRestarts_;

  /** max number of iterations */
  int maxIterations_;

  /** convergence tolerance */
  double tol_;
};

/**
 *  @class FieldImplicitDirectEulerIntegrator
 *  @brief implicit Euler method with direct solve
 */
class FieldImplicitDirectEulerIntegrator : public FieldEulerIntegrator {

 public:

  /** Constructor */
  FieldImplicitDirectEulerIntegrator(
    const FieldName fieldName,
    const PhysicsModel * physicsModel,
    /*const*/ FE_Engine * feEngine,
    /*const*/ ATC_Coupling * atc,
    const Array2D< bool > & rhsMask, // copy
    const double alpha = 0.5 // default to trap/midpt
  );

  /** Destructor */
  virtual ~FieldImplicitDirectEulerIntegrator();

  /** initalize - init the matrices and inverses */
  void initialize(const double dt, const double time,
    FIELDS & fields);

  /** update */
  void update(const double dt, const double time,
    FIELDS & fields,  FIELDS & rhs);

 protected:
  /** euler update factor */
  double alpha_;

  /** matrices */
  SPAR_MAT  _M_;
  SPAR_MAT  _K_;
  SPAR_MAT  _lhsMK_;
  SPAR_MAT  _rhsMK_;
  class LinearSolver * solver_;
};

} // namespace ATC

#endif
