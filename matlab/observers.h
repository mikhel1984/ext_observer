#ifndef OBSERVERS_DYN_H
#define OBSERVERS_DYN_H

#define JOINT_NO 2

#define ERR_WRONG_INDEX -1
#define ERR_NO_SLOTS    -2
#define ERR_WRONG_TYPE  -3

#define ADD_NEW -1
// If index -1, return id of the new observer, else update existant

#ifdef __cplusplus
extern "C" {
#endif


int configMomentumObserver(int index, double k[JOINT_NO]);

int configDisturbanceObserver(int index, double sigma, double xeta, double beta);

int configSlidingModeObserver(int index, double T1[JOINT_NO], double S1[JOINT_NO], double T2[JOINT_NO], double S2[JOINT_NO]);

int configDistKalmanObserver(int index, double S[JOINT_NO*JOINT_NO], double H[JOINT_NO*JOINT_NO], double Q[4*JOINT_NO*JOINT_NO], double R[JOINT_NO*JOINT_NO]);

int configDistKalmanObserverExp(int index, double S[JOINT_NO*JOINT_NO], double H[JOINT_NO*JOINT_NO], double Q[4*JOINT_NO*JOINT_NO], double R[JOINT_NO*JOINT_NO]);

int configFilterDynObserver(int index, double cutOff, double dt);

int configFilterRangeObserver(int ind, double cutOff, double dt, double k);

// Get external torque estimation
int getExternalTorque(int index, double ext[JOINT_NO], double q[JOINT_NO], double qd[JOINT_NO], double tau[JOINT_NO], double dt);

// reset observer state
void reset(int index);

// clear memory
void freeAll();

// get expected joint torques
void getRobotTorque(double tau[JOINT_NO], double q[JOINT_NO], double qd[JOINT_NO], double q2d[JOINT_NO]);

#ifdef __cplusplus
}
#endif


#endif // OBSERVERS_DYN_H
