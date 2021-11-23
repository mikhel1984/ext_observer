# ext_observer

External torque observers for a robot manipulator. 

**lib** contains C++ implementations for such observers as:
- momentum observer
- disturbance observer
- sliding mode observer 
- Kalman filter based observer
- dynamic filter based observer

**matlab** is the Matlab wrapper code for the observers.

**tests** contains some tests.

Integration with **ROS** can be found in the "ros_examples" branch.

The library was developed for the research 
```
Mamedov, S., Mikhel, S. (2020). Practical aspects of model-based collision detection. 
Frontiers in Robotics and AI, 7, 162.
```
