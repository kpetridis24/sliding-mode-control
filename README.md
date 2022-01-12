# sliding-mode-control
Implementation of Sliding Mode Control and testing on both regulation and tracking problems

This is a MATLAB implementation of the [Sliding Mode Control](https://en.wikipedia.org/wiki/Sliding_mode_control) method for controlling dynamic systems,
which include uncertainties. This means that we know the ODE of the control system, but we don't acquire any knowledge about the exact values of the system's parameters.
For the exhibition of the method, i simulated the dynamic system of a robotic arm, with 2 freedom degrees, and designed a controller following the Sliding Mode Control 
method. The controller is tested on both a regulation and a tracking problem. Regulation refers to us, dictating that the system's output is going to converge into a specific point and remain there. Tracking is the problem where, we desire to make the system track a specific orbit with high accuracy, after some time t.
