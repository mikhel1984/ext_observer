/**
 * @file external_observer.h
 *
 * @brief Abstract classes for a robot and external torque observer.
 *
 * Define robot interfaces for working with dynamic matrices and recursice Newton-Euler algorithm
 * and corresponding observers. 
 */
#ifndef EXTERNAL_OBSERVER_H
#define EXTERNAL_OBSERVER_H

#include <eigen3/Eigen/Geometry>

/** @brief Gravity constant. */
#define GRAVITY 9.81

/**
 * @brief Simplified object call.
 */
typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

/**
 * @brief Common elements of a robot interface.
 */
class RobotDynamicsBase {
public:
  /** 
   * @brief Friction model.
   * 
   * Find expected friction using velocity.
   * @param qd vector of joint velocities
   * @return torque loss due to friction.
   */
  virtual Vector getFriction(Vector& qd) = 0;
  /** 
   * @brief Number of movable joints in robot.
   * @return Joint number.
   */
  virtual int jointNo() = 0;

}; // RobotDynamicsBase

/** 
 * @brief Base class when M, C and G are known.
 * 
 * Base class for a robot if all dynamical matrices could be found explicitly.
 */ 
class RobotDynamics : public RobotDynamicsBase {
public:
  /** 
   * @brief Get inertia matrix.
   * @param q joint angle vector.
   * @return matrix M.
   */
  virtual Matrix getM(Vector& q) = 0;
  /** 
   * @brief Get Coriolis/centrifugal matrix.
   * @param q joint angle vector.
   * @param qd joint angle velocity vector.
   * @return matrix C.
   */
  virtual Matrix getC(Vector& q, Vector& qd) = 0;
  /** 
   * @brief Get gravity terms.
   * @param q joint angle vector. 
   * @return vector G.
   */
  virtual Vector getG(Vector& q) = 0;
  
}; // RobotDynamics

/**
 * @brief Base class when RNEA is defined.
 * 
 * Base class for a robot if only RNEA could be applied. 
 */
class RobotDynamicsRnea : public RobotDynamicsBase {
public:
  /**
   * @brief Object constructor.
   */
  RobotDynamicsRnea()
  : RobotDynamicsBase()
  , _qext(Vector(1))
  , _p0(Vector(1))
  , _zero(Vector(1))
  , _sum(Vector(1))
  , _M(Matrix(1,1))
  , delta(1E-7) 
  {   }
  /** 
   * @brief Call recursibe Newton-Euler algorithm.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @param q2d joint acceleraiton vector.
   * @param g gravity constant.
   * @return torque value.
   */
  virtual Vector rnea(Vector& q, Vector& qd, Vector& q2d, double g=0) = 0;
  /** 
   * @brief Angle step for derivative estimation.
   * 
   * Joint angle step for numerical evaluation of dM/dt. 
   * @param d desired value (default is ~ 1E-7).
   */
  void setDelta(double d) { delta = d; }
  /**
   * @brief Estimation of product C^T * qd.
   *
   * Use numerical differentiation technique. 
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @return product value estiamtion.
   */
  Vector tranCqd(Vector& q, Vector& qd);
  /**
   * @brief Build matrix of inertia.
   * @param q joint angle vector.
   * @return matrix M.
   */
  Matrix getM(Vector& q);
  
protected:
  /** @brief Temporary variables, avoid memory reallocation. */
  Vector _qext, _p0, _zero, _sum; 
  Matrix _M; 
  /** @bried Differentiation step. */
  double delta;
  
}; // RobotDynamicsRnea

/** 
 * @brief External torque observer.
 *
 * Work with matrices M, C ang G.
 */
class ExternalObserver {
public:
  /**
   * @brief Object constructor. 
   * @param rd pointer to the robot interface.
   */
  ExternalObserver(RobotDynamics *rd)
  : dyn(rd) 
  , isRun(false)
  { jointNo = rd->jointNo(); }
  /**
   * @brief Estimate external torques.
   * 
   * The main method to use.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @param tau joint torque vector.
   * @param dt time step from the previous call.
   * @return external torque vector.
   */
  virtual Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt) = 0;
  /** 
   * @brief Reset observer state.
   */ 
  void reset() { isRun = false; }

protected:
  RobotDynamics *dyn;    /**< Pointer to the robot interface. */
  int jointNo;           /**< Number of joints. */
  bool isRun;            /**< First call check. */

}; // ExternalObserver

/**
 * @brief External torque observer.
 * 
 * Work with RNEA algorithm.
 */
class ExternalObserverRnea {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot interface.
   */
  ExternalObserverRnea(RobotDynamicsRnea *rd) 
  : dyn(rd)
  , isRun(false)
  { jointNo = rd->jointNo(); }
  /**
   * @brief Estimate external torques.
   * 
   * The main method to use.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @param tau joint torque vector.
   * @param dt time step from the previous call.
   * @return external torque vector.
   */
  virtual Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt) = 0;
  /** 
   * @brief Reset observer state.
   */ 
  void reset() { isRun = false; }  

protected:
  RobotDynamicsRnea *dyn; /**< Pointer to the robot interface. */
  int jointNo;            /**< Number of joints. */
  bool isRun;             /**< First call check. */

}; // ExternalObserverRnea

// Find C^T * qd
Vector RobotDynamicsRnea::tranCqd(Vector& q, Vector& qd)
{
  // TODO: call it once
  int N = jointNo();
  _zero.resize(N); _zero.setZero();
  _p0 = rnea(q,_zero,qd);  // M * qd
  _sum.setZero();

  for(int i = 0; i < N; i++) {
    _qext = q;
    _qext(i) += delta; 
    _sum += (rnea(_qext,_zero,qd) - _p0) * (qd(i) / delta);
  }
  _sum -= rnea(q,qd,_zero); // M'*qd - C*qd

  return _sum;
}

// Create M matrix from sequence of RNEA calls
Matrix RobotDynamicsRnea::getM(Vector& q)
{
  // TODO: call it once
  int N = jointNo();
  _zero.resize(N); _zero.setZero();
  _M.resize(N,N); 
  _qext.resize(N);
  
  for(int i = 0; i < N; i++) {
    _qext.setZero();  // use _qext to choose column
    _qext(i) = 1;   
    _p0 = rnea(q,_zero,_qext); // use _p0 to save matrix column
    for(int j = 0; j < N; j++)
      _M(j,i) = _p0(j);
  }
  
  return _M;
}

#endif // EXTERNAL_OBSERVER_H

