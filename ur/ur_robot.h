#ifndef UR_ROBOT_H
#define UR_ROBOT_H

#include "../lib/external_observer.h"

#define NJ 6
#define NPARAM 60
#define NFRIC 18
#define NM 36

class UrRobot : public RobotDynamics {
public:
  UrRobot();
  ~UrRobot();

  int jointNo() { return NJ; }

  Matrix getM(Vector& q);
  
  Matrix getC(Vector& q, Vector& qd);
  
  Vector getG(Vector& q);
  
  Vector getFriction(Vector& qd);


private:
  Matrix M, C;
  Vector G, F;
  
  double mArr[NM], gArr[NJ], cArr[NM], fArr[NJ];
  double qArr[NJ], qdArr[NJ];
  double *param, *fric, *dyn;
};

#endif // UR_ROBOT_H
