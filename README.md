UR-5e Forward and Inverse Kinematics Analysis
📌 Project Overview

This project implements a kinematic analysis of the Universal Robots UR-5e, a 6-DOF industrial manipulator. It derives the Forward Kinematics (FK) using standard Denavit-Hartenberg (DH) parameters and solves the Inverse Kinematics (IK) numerically using the Newton-Raphson method.

The core objective is to simulate the robot following a constrained elliptical trajectory in task space, visualized via a 3D "stick figure" animation in MATLAB.
🚀 Key Features

    DH Parameter Derivation: Calculation of the homogeneous transformation matrices from the base to the end-effector.

    Numerical Inverse Kinematics: Implementation of the Newton-Raphson iterative algorithm to resolve joint angles for target poses.

    Trajectory Generation: Parametric equation setup for an elliptical path in 3D space.

    MATLAB Visualization: A frame-by-frame 3D animation of the robot links and joints tracking the path.

⚙️ Methodology
1. Forward Kinematics

The forward kinematics are computed by assigning coordinate frames to each of the 6 links. The overall transformation matrix is derived as:
0T6​=0T1​⋅1T2​⋅2T3​⋅3T4​⋅4T5​⋅5T6​
2. Inverse Kinematics (Newton-Raphson)

Since the UR-5e has a complex workspace, a numerical approach is used to find the joint configuration q that minimizes the error between the current end-effector position and the target trajectory point.

The update rule used is:
qk+1​=qk​+J†(qk​)⋅(xtarget​−xcurrent​)

Where J† represents the pseudo-inverse of the Jacobian.
📂 Repository Structure
Plaintext
