// Copyright 2020-2024 Stanislav Mikhel

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
const double GRAVITY = 9.81;

/**
 * @brief Simplified object call.
 */
typedef Eigen::MatrixXd MatrixJ;
typedef Eigen::VectorXd VectorJ;

/**
 * @brief Common elements of a robot interface.
 */
class RobotDynamicsBase {
public:
  virtual ~RobotDynamicsBase() = default;

  /**
   * @brief Friction model.
   *
   * Find expected friction using velocity.
   * @param qd vector of joint velocities
   * @return torque loss due to friction.
   */
  virtual VectorJ getFriction(VectorJ& qd) = 0;

  /**
   * @brief Number of movable joints in robot.
   * @return Joint number.
   */
  virtual int jointNo() = 0;
};  // RobotDynamicsBase


/**
 * @brief Base class when M, C and G are known.
 *
 * Base class for a robot if all dynamical matrices could be found explicitly.
 */
class RobotDynamics : public RobotDynamicsBase {
public:
  virtual ~RobotDynamics() = default;

  /**
   * @brief Get inertia matrix.
   * @param q joint angle vector.
   * @return matrix M.
   */
  virtual MatrixJ getM(VectorJ& q) = 0;

  /**
   * @brief Get Coriolis/centrifugal matrix.
   * @param q joint angle vector.
   * @param qd joint angle velocity vector.
   * @return matrix C.
   */

  virtual MatrixJ getC(VectorJ& q, VectorJ& qd) = 0;
  /**
   * @brief Get gravity terms.
   * @param q joint angle vector.
   * @return vector G.
   */
  virtual VectorJ getG(VectorJ& q) = 0;
};  // RobotDynamics


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
  RobotDynamicsRnea();

  virtual ~RobotDynamicsRnea() = default;

  /**
   * @brief Call recursibe Newton-Euler algorithm.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @param q2d joint acceleraiton vector.
   * @param g gravity constant.
   * @return torque value.
   */
  virtual VectorJ rnea(VectorJ& q, VectorJ& qd, VectorJ& q2d, double g = 0) = 0;

  /**
   * @brief Angle step for derivative estimation.
   *
   * Joint angle step for numerical evaluation of dM/dt.
   * @param d desired value.
   */
  inline void setDelta(double d) noexcept { delta = d; }

  /**
   * @brief Estimation of product C^T * qd.
   *
   * Use numerical differentiation technique.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @return product value estiamtion.
   */
  VectorJ tranCqd(VectorJ& q, VectorJ& qd);

  /**
   * @brief Build matrix of inertia.
   * @param q joint angle vector.
   * @return matrix M.
   */
  MatrixJ getM(VectorJ& q);

protected:
  /** @brief Temporary variables, avoid memory reallocation. */
  VectorJ _qext, _p0, _zero, _sum;
  MatrixJ _M;
  /** @bried Differentiation step. */
  double delta = 1E-7;
};  // RobotDynamicsRnea


class ExternalObserverBase {
public:
  /**
   * @brief Object constructor.
   */
  ExternalObserverBase()
  : jointNo(0)
  , objType(-1)
  , isRun(false) {}

  virtual ~ExternalObserverBase() = default;

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
  virtual VectorJ getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt) = 0;

  /**
   * @brief Reset observer state.
   */
  inline void reset() noexcept { isRun = false; }

  /**
   * @brief Get observer type.
   */
  inline int type() const noexcept { return objType; }

protected:
  int jointNo;           /**< Number of joints. */
  int objType;           /**< Observer type. */
  bool isRun;            /**< First call check. */
};  // ExternalObserverBase


/**
 * @brief External torque observer.
 *
 * Work with matrices M, C ang G.
 */
class ExternalObserver : public ExternalObserverBase {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot interface.
   */
  ExternalObserver(RobotDynamics *rd, int type)
  : ExternalObserverBase()
  , dyn(rd)
  {
    jointNo = rd->jointNo();
    objType = type;
  }

  virtual ~ExternalObserver() = default;

protected:
  RobotDynamics *dyn;    /**< Pointer to the robot interface. */
};  // ExternalObserver


/**
 * @brief External torque observer.
 *
 * Work with RNEA algorithm.
 */
class ExternalObserverRnea : public ExternalObserverBase {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot interface.
   */
  ExternalObserverRnea(RobotDynamicsRnea *rd, int type)
  : ExternalObserverBase()
  , dyn(rd)
  {
    jointNo = rd->jointNo();
    objType = type;
  }

  virtual ~ExternalObserverRnea() = default;

protected:
  RobotDynamicsRnea *dyn; /**< Pointer to the robot interface. */
};  // ExternalObserverRnea

#endif  // EXTERNAL_OBSERVER_H

