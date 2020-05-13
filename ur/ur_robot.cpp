
#include "ur_robot.h"

#include "C_mtrx_fcn/C_mtrx_fcn.h"
#include "C_mtrx_fcn/C_mtrx_fcn_initialize.h"
#include "C_mtrx_fcn/C_mtrx_fcn_terminate.h"

#include "F_vctr_fcn/F_vctr_fcn.h"
#include "F_vctr_fcn/F_vctr_fcn_initialize.h"
#include "F_vctr_fcn/F_vctr_fcn_terminate.h"

#include "G_vctr_fcn/G_vctr_fcn.h"
#include "G_vctr_fcn/G_vctr_fcn_initialize.h"
#include "G_vctr_fcn/G_vctr_fcn_terminate.h"

#include "M_mtrx_fcn/M_mtrx_fcn.h"
#include "M_mtrx_fcn/M_mtrx_fcn_initialize.h"
#include "M_mtrx_fcn/M_mtrx_fcn_terminate.h"

double drvsLst[NJ] = {0.0000,3.7434,0.0000,0.0708,0.2342,0.4365}; 

double fricLst[NFRIC] = {21.256003,12.548391,0.196127,20.220907,
13.265644,-0.753542,10.377742,4.994313,0.197915,3.575050,1.997854,
0.054820,2.490901,2.690460,-0.0096109,3.025553,2.297854,0.042528};

double rgdLst[NPARAM] = {16.017661,-2.0048e-24,-4.2449e-22,16.017661,
-4.2449e-22,9.4372e-07,-6.5344e-22,-6.5344e-22,3.50506e-15,
4.827358,3.932521,0.852873,0.151907,3.720556,-0.594538,1.903157,
0.049952,-0.275933,3.263268,7.925888,2.224532,-0.536564,-0.176607,
1.388985,-0.216689,1.388684,0.023948,0.242903,0.893879,2.480855,
1.029737,0.232810,-0.012259,0.075126,0.074354,1.065273,-0.000583,
-0.253394,0.191763,2.155997,0.028615,-0.033833,-0.025661,0.108051,
-0.013845,0.111089,0.007875,0.093425,0.114659,2.155994,0.028384,
-0.002076,-0.013982,0.035360,-0.004008,0.008145,0.002666,-0.001080,
-0.015801,0.222170};

double drvsDif[NJ] = {0, 0.0889, 0.0789, 0.0618, 0.0651, 0.0567};
double fricDif[NFRIC] = {0.0686, 0.0410, 0.0222, 0.0949, 0.0451, 0.0759, 0.0985,
0.0378, 0.0512, 0.0865, 0.0380, 0.0240, 0.0862, 0.0382,
0.0232, 0.0744, 0.0421, 0.0233};
double rgdDif[NPARAM] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0502, 0.0411, 0, 
0.0542, 0.0766, 0.0085, 0, 0.0052, 0, 0.0765, 0.0454,
0.0355, 0.0604, 0.0410, 0, 0.0072, 0, 0.0054, 0, 0.0758,
0.0318, 0.0313, 0, 0.0336, 0.0719, 0.0040, 0, 0, 0, 
0.0621, 0.0171, 0.0245, 0.0425, 0.0219, 0, 0.0041, 0.0030,
0.0033, 0, 0.0423, 0.0159, 0.0144, 0.0274, 0.0146, 0.0429,
0.0020, 0, 0.0026, 0};


UrRobot::UrRobot()
       : RobotDynamics()
       , M(Matrix(NJ,NJ))       
       , C(Matrix(NJ,NJ))
       , G(Vector(NJ))
       , F(Vector(NJ))
       , param(rgdLst)
       , fric(fricLst)
       , dyn(drvsLst)
{
  C_mtrx_fcn_initialize();
  F_vctr_fcn_initialize();
  G_vctr_fcn_initialize();
  M_mtrx_fcn_initialize();
  // parameters 
        
  for(int i = 0; i < NJ; i++) drvsDif[i] += drvsLst[i];
  for(int i = 0; i < NFRIC; i++) fricDif[i] += fricLst[i];
  for(int i = 0; i < NPARAM; i++) rgdDif[i] += rgdLst[i];   
}

UrRobot::~UrRobot()
{
  C_mtrx_fcn_terminate();
  F_vctr_fcn_terminate();
  G_vctr_fcn_terminate();
  M_mtrx_fcn_terminate();
}

Matrix UrRobot::getM(Vector& q)
{
  for(int i = 0; i < NJ; i++) qArr[i] = q(i);
  // evaluate
  M_mtrx_fcn(qArr,param,mArr); 
  // save
  for(int i = 0; i < NJ; i++){
    for(int j = 0; j < NJ; j++) 
      M(i,j) = mArr[NJ*j+i];
  }
  for(int i = 0; i < NJ; i++) 
    M(i,i) += dyn[i];
  
  return M;
}

Matrix UrRobot::varM(Vector& q)
{
  for(int i = 0; i < NJ; i++) qArr[i] = q(i);
  // M + dM
  M_mtrx_fcn(qArr,rgdDif,mArr); 
  // save
  for(int i = 0; i < NJ; i++){
    for(int j = 0; j < NJ; j++) 
      M(i,j) = mArr[NJ*j+i];
  }
  for(int i = 0; i < NJ; i++) 
    M(i,i) += drvsDif[i];
  
  // M 
  M_mtrx_fcn(qArr,param,mArr); 
  
  for(int i = 0; i < NJ; i++){
    for(int j = 0; j < NJ; j++) 
      M(i,j) -= mArr[NJ*j+i];
  }
  for(int i = 0; i < NJ; i++) 
    M(i,i) -= dyn[i];
  
  return M;
}

Vector UrRobot::getG(Vector& q)
{
  for(int i = 0; i < NJ; i++) qArr[i] = q(i);
  // evaluate
  G_vctr_fcn(qArr,param,gArr);
  // save 
  for(int i = 0; i < NJ; i++) 
    G(i) = gArr[i];
  return G;
}

Vector UrRobot::varG(Vector& q)
{
  for(int i = 0; i < NJ; i++) qArr[i] = q(i);
  // G + dG
  G_vctr_fcn(qArr,rgdDif,gArr);
  // save 
  for(int i = 0; i < NJ; i++) 
    G(i) = gArr[i];
        
  // G
  G_vctr_fcn(qArr,param,gArr);
  // save 
  for(int i = 0; i < NJ; i++) 
    G(i) -= gArr[i];
        
  return G;
}

Matrix UrRobot::getC(Vector& q, Vector& qd)
{
  for(int i = 0; i < NJ; i++) {
    qArr[i] = q(i);
    qdArr[i] = qd(i);
  }
  // evaluate
  C_mtrx_fcn(qArr,qdArr,param,cArr);
  // save
  for(int i = 0; i < NJ; i++){
    for(int j = 0; j < NJ; j++) 
      C(i,j) = cArr[NJ*j+i];
  }
  return C; 
}

Matrix UrRobot::varC(Vector& q, Vector& qd)
{
  for(int i = 0; i < NJ; i++) {
    qArr[i] = q(i);
    qdArr[i] = qd(i);
  }
  // C + dC
  C_mtrx_fcn(qArr,qdArr,rgdDif,cArr);
  // save
  for(int i = 0; i < NJ; i++){
    for(int j = 0; j < NJ; j++) 
      C(i,j) = cArr[NJ*j+i];
  }
  
  // C
  C_mtrx_fcn(qArr,qdArr,param,cArr);
  // save
  for(int i = 0; i < NJ; i++){
    for(int j = 0; j < NJ; j++) 
      C(i,j) -= cArr[NJ*j+i];
  }
  
  return C; 
}

Vector UrRobot::getFriction(Vector& qd)
{
  for(int i = 0; i < NJ; i++) qdArr[i] = qd(i);
  // evaluate
  F_vctr_fcn(qdArr,fric,fArr);
  // copy
  for(int i = 0; i < NJ; i++) 
    F(i) = fArr[i];
  return F;
}

Vector UrRobot::varFriction(Vector& qd)
{
  for(int i = 0; i < NJ; i++) qdArr[i] = qd(i);
  // evaluate
  F_vctr_fcn(qdArr,fricDif,fArr);
  // copy
  for(int i = 0; i < NJ; i++) 
    F(i) = fArr[i];
        
  F_vctr_fcn(qdArr,fric,fArr);
  // copy
  for(int i = 0; i < NJ; i++) 
    F(i) -= fArr[i];
        
  return F;
}
