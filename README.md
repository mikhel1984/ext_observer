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

The library was developed for the [research](https://www.frontiersin.org/articles/10.3389/frobt.2020.571574/full?utm_campaign=Artificial%2BIntelligence%2BWeekly&utm_medium=web&utm_source=Artificial_Intelligence_Weekly_189) 
```
@article{mamedov2020practical,
  title={Practical aspects of model-based collision detection},
  author={Mamedov, Shamil and Mikhel, Stanislav},
  journal={Frontiers in Robotics and AI},
  volume={7},
  pages={571574},
  year={2020},
  publisher={Frontiers Media SA}
}
```
