# Implementation of the torque observers 

The main method for each observer is 
```cpp
Vector getExternalTorque(Vector& q, Vector& dq, Vector& tau, double dt)
```
where _q_ is the vector of joint angles, _dq_ is the vector of joint velocities, _tau_ is the vector of measured joint torques, _dt_ is the time step. 

The robot dynamics can be defined in 2 ways. If all the matrices (_M_, _C_ and _G_) are known exactly, use default implementation. When dynamics is defined in form of RNEA algorithm (i.e. as a function _f(q,dq,ddq)_) use implementations with "Rnea" in the end.
