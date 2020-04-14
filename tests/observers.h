#ifndef OBSERVERS_DYN_H
#define OBSERVERS_DYN_H 

#define ERR_WRONG_INDEX -1
#define ERR_NO_SLOTS    -2
#define ERR_WRONG_TYPE  -3

// If index -1, return id of the new observer, else update existant

int configMomentumObserver(int index, double k[2]);

int configDisturbanceObserver(int index, double sigma, double xeta, double beta);

int configSlidingModeObserver(int index, double T1[2], double S1[2], double T2[2], double S2[2]);

int configDistKalmanObserver(int index, double S[2*2], double H[2*2], double Q[4*4], double R[2*2]);

int configFilterDynObserver(int index, double cutOff, double dt);

// Get external torque estimation
int getExternalTorque(int index, double ext[2], double q[2], double qd[2], double tau[2], double dt); 

// reset observer state
void reset(int index);

// clear memory
void freeAll();


#endif // OBSERVERS_DYN_H
