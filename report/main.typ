#import "@preview/subpar:0.2.2"
#import "@preview/physica:0.9.8": *
#import "@preview/clean-cnam-template:1.5.0": *
#import "@preview/bubble:0.2.2": *

#show: super-T-as-transpose 

#set document(title: "")
#show math.mat: math.display
#show math.frac: math.display

#show: bubble.with(
  title: "Energy-Aware Path Choice with Computed-Torque Tracking",
  subtitle: "Robotics Project",
  author: "Ahmad Abu Zainab - 6320",
  affiliation: "Lebanese University  - Faculty of Engineering III",
  date: "2025-2026",
  // year: "2025-2026",
  class: "Robotics",
  // other: ("Made with Typst", "https://typst.com"),
  main-color: "2e63a9", //set the main color
  logo: image("figs/ulfg.jpg", width: 30%),
) 

#outline(title: "Table of Contents")
#pagebreak()

= Introduction

The goal of this project is to track different end-effector paths to compute torque (via inverse dynamics), tracking error and total energy consumed at each path. For this project, we will be using a simplified model of the SCARA robot configuration and resolved rate control with PID for control as resolved rate control doesn't involve forward dynamics as other control methods do. #cite(<whitney2007resolved>)


= Theory and Implementation

For this project we will use a simplified version of the SCARA robot model in which we make a few assumptions to simplify our work
+ No joint limits (This doesn't really change much as we can simply take a path that doesnt require us to exceed usual SCARA robot joint limits)
+ A movement of $Delta q$ of the robot configuration is always possible and unimpeded (If we command the robot to move to a certain position in the $q$-space it will always perform that movement)


== Kinematics

=== DH Parameters

First we define the robot parameters as follows

#subpar.grid(
  figure(image("figs/structure.png"), caption: [The structure of the SCARA robot]),
  figure(image("figs/coordinate.png"), caption: [The coordinate systems of the SCARA robot]),
  columns: (1fr, 1fr),
  caption: [Configuration of the SCARA Robot #cite(<feng2016kinematic>)],
  label: <config>,
)


#figure(
  table(
    columns: 5,
    inset: 7pt,
    stroke: none,
    gutter: 2pt,
    table.hline(),
    table.header([Joints], [$d_i$], [$a_(i-1)$], [$alpha_(i-1)$], [$theta_i$]),
    table.hline(stroke: (dash:"dashed")),
    [1], [$1$], [$0.5$], [$0$], [$theta_1$],
    [2], [$0$], [$0.5$], [$0$], [$theta_2$],
    [3], [$-d_3$], [$0$], [$0$], [$0$],
    [4], [$0$], [$0$], [$0$], [$theta_4$],
    table.hline(),
  ),
  caption: [DH Parameters of the SCARA Robot]
)

We apply the non-proximal DH parameter transformation matrix. #cite(<denavit1955kinematic>)

$
attach(T, tl: i-1, bl: i) = mat(
  delim: "[",
  cos theta_i, -sin theta_i cos alpha_i, sin theta_i sin alpha_i, a_i cos theta_i;
  sin theta_i, cos theta_i cos alpha_i, -cos theta_i sin alpha_i, a_i sin theta_i;
  0, sin alpha_i , cos alpha_i, d_i;
  0, 0, 0, 1  
)
$

Here we can lock the orientation of the end-effector ($theta_4=0$) as it's not really relevant to our study, so our $vb(q)$ vector becomes.

$
vb(q) = mat(delim:"[", theta_1, theta_2, d_3)^T
$

Finally we obtain the final transformation matrix:

$
attach(T, tl: 3, bl: 0) = mat(
  delim: "[",
  c_(124), -s_(124), 0, 0.5 c_(1) + 0.5c_(12);
  s_(124), c_(124), 0, 0.5 s_(1) + 0.5s_(12);
  0, 0 , 1, 1-d_3;
  0, 0, 0, 1  
)
$

=== Forward Kinematics

Using the transformation matrix we can get the position of the end-effector $vb(x)$

$
vb(x) = f(vb(q)) = mat(delim:"[", 0.5 cos theta_1 + 0.5cos theta_1+theta_2 ; 0.5 sin theta_1 + 0.5sin theta_1+theta_2 ; 1-d_3)
$


=== Inverse Kinematics

While inverse kinematics are not needed in this study we can easily deduce it from the forward kinematics

$
vb(q) = mat(delim:"[", op("atan2")(x,y) - op("atan2")(1/2+gamma/2,delta/2); plus.minus op("acos")(2(x^2+y^2)-1) ; 1-z) "where" gamma = 2x^2 + 2y^2 - 1 "and" delta = sqrt(1-gamma^2)
$

== Jacobian

We can easily obtain the Jacobian and inverse Jacobian using direct diffrentiation

$
vb(J) = mat(delim:"[",
pdv(f_1,q_1),pdv(f_1,q_2) ,..., pdv(f_1,q_n);
dots.v,dots.v,dots.down,dots.v,;
pdv(f_m,q_1),pdv(f_m,q_2) ,..., pdv(f_m,q_n);
) = mat(delim:"[",-0.5sin(θ_1) - 0.5sin(θ_1 + θ_2), -0.5sin(θ_1 + θ_2),0; 0.5cos(θ_1) + 0.5cos(θ_1 + θ_2), 0.5cos(θ_1 + θ_2),0;0,0,-1)
\
vb(J)^(-1) = mat(delim: "[", 
frac(2 cos (θ_(1) + θ_(2)), sin (θ_(2))), frac(2 sin (θ_(1) + θ_(2)), sin (θ_(2))), 0 ;
- frac(2 cos (θ_(1)) + 2 cos (θ_(1) + θ_(2)), sin (θ_(2))), - frac(2 sin (θ_(1)) + 2 sin (θ_(1) + θ_(2)), sin (θ_(2))), 0;
0, 0, - 1) 
$

== Selected Trajectory

For this project we selected 3 (reachable) points arbitrarily 
$
vb(p_1) = f(-pi/6, -pi/2, 0) = mat(
  delim:"[",
  0.1830127;
  -0.6830127;
  1
) 
quad 
vb(p)_2 = f(pi/6, pi/2, 1) = mat(
  delim:"[",
  0.1830127;
  0.6830127;
  0
)\
vb(p)_i = mat(
  delim:"[",
  0.5;
  0;
  0.75
)
$

Then we define 3 paths($t in [0,1]$): a straight line from $vb(p)_1$ to $vb(p)_2$, a straight line from $vb(p)_1$ to an intermediate point $vb(p)_i$ then to $vb(p)_2$, and a cubic spline from $vb(p)_1$ to $vb(p)_2$ passing through $vb(p)_i$ ($vb(P)'_3 "at" vb(p_i) = f "otherwise" 0$).

$
vb(P)_1(t) = (vb(p)_2 - vb(p)_1)t + vb(p)_1 quad vb(P)_2(t) = cases(
  2(vb(p)_i - vb(p)_1)t + vb(p)_1 &"if" t<= 0.5,
  2(vb(p)_2 - vb(p)_i)(t-0.5) + vb(p)_i &"if" t> 0.5
)\
vb(P)_3(t) = cases(
  (2vb(p)_1 - 2vb(p)_i+f)2^3t^3 + (3vb(p)_i - 3vb(p)_1 - f)2^2t^2 + vb(p)_1 &"if" t<= 0.5,
  (2vb(p)_i - 2vb(p)_2+f)2^3(t-0.5)^3 + (3vb(p)_2 - 3vb(p)_i - 2f)2^2(t-0.5)^2 + 
  2f(t-0.5) + vb(p)_i &"if" t> 0.5
)
$

== Dynamics

While the dynamics for the SCARA robot are readily available they are often large and cumbersome to implement, so for this project we can derive our own from a simplified model of the SCARA robot by assuming the mass of the robot is concentrated at the joints. So for joint 1 the mass is $m_1=5"kg"$, joint 2 $m_2=4"kg"$, and joint 3 $m_3=2"kg"$.

We can the compute the K.E. and P.E. for each joint

$
T = sum_(i=1)^m (1/2 m_i display(sum_(j=1)^n (dv(q_j,t))^2 ))\
U = sum_(i=1)^m g dot z_i dot m_i\
L = T-U
$

Finally we plug in $L$ in the Euler–Lagrange equations #cite(<ortega1998euler>)

$
phi_i = dv(,t) pdv(L, dot(q)_i) - pdv(L, q_i)
$

Noting that $phi_i$ is the generalized force where

$
phi_i = cases(
  tau_i & "if the joint is revolute",
  f_i & "if the joint is prismatic"
)
$

Computing the resultant torques is tedious, so we can use a package like `sympy` to compute the torques for us #cite(<sympy>). Finally we obtain the following dynamic equation

$
vb(phi) = mat(
  delim: "[",
  - 3 sin (θ_(2)) dot(theta)_(1) dot(theta)_(2) - 1.5 sin (θ_(2))dot(theta)_(2)^(2) + 3 cos(θ_(2)) dot.double(theta)_(1) + 1.5 cos (θ_(2)) dot.double(theta)_(2) + 4.25 dot.double(theta)_(1) + 1.5 dot.double(theta)_(2);
    1.5 sin(θ_(2)) dot(theta)_(1)^(2) + 1.5 cos(θ_(2)(t)) dot.double(theta)_(1) + 1.5 dot.double(theta)_(1) + 1.5 dot.double(theta)_(2);
    2  dot.double(d_3) - 19.6
)
$

== Control

As mentioned in the introduction, we will use  resolved rate control with PID control. The PID parameters are as the following

$
K_p = 10.0 quad K_i = 0.5 quad K_d = 0.05
$

#figure(
  image("figs/control.png", width: 85%),
  caption: [Resolved Rate Control Scheme #cite(<whitney2007resolved>)]
)

== Implementation 

The simulation was implemented in Python using libraries like `numpy` for computation and `matplotlib` for plots #cite(<mat>) #cite(<numpy>). `vpython` was used for the 3D visualization. The code is bundled in Python notebooks. The project contains 2 Python notebooks `Derivation.ipynb` and `Simulation.ipynb`.
- `Derivation.ipynb` contains the code used to derive the transformation matrix, Jacobian, dynamics, and the position of every joint along the robot for the 3D visualization.
- `Simulation.ipynb` contains the code for the simulation and visualization and it's functions are a direct implementation of the ones derived in `Derivation.ipynb`.

To select the path in the code set the `path` variables to either `path1`, `path2`, or `path3`. 

=== Visualization

Upon selecting the path and running all the cells below the cell that selects the path, you can see that the last cell contains a `vpython` window showing a robot visualization which can be rotated around by dragging it using the Right Mouse Button. Grey cylinders are fixed parts/support of the robot, red spheres are joints, blue cylinders are links and the green sphere is the end-effector. The purple path is the desired path for the robot.

#figure(
  image("figs/robot-vpython.png", width: 30%),
  caption: [SCARA Robot Visualization in `vpython`]
)



= Results and Analysis

Running the simulation for all 3 paths we obtain the following measurements for energy

#figure(
  table(
    columns: 2,
    stroke: none,
    align: center,
    table.hline(),
    table.header([Path Function], [Energy Consumed]),
    table.hline(stroke: (dash:"dashed")),
    [$vb(P)_1(t)$], [$19.1803 "J"$],
    [$vb(P)_2(t)$], [$18.8058 "J"$],
    [$vb(P)_3(t)$], [$37.2540 "J"$],
    table.hline(),
  ),
  caption: [Resultant Energy for Each Path]
)

We notice that path 2 consumes the least energy despite not being the shortest path.

#subpar.grid(
  figure(image("figs/J1.png")),
  figure(image("figs/J2.png")),
  grid.cell(figure(image("figs/J3.png", width: 50%)), colspan: 2),
  columns: (1fr, 1fr),
  caption: [Tracking Error for Each Path],
  label: <joints>,
)

We can see from each plot what kind of discontinuity each path has. Path 1 is a straight line so it doesn't have any discontinuity and is $C^infinity$ continuos, path 2, on the other hand, clearly has a sharp edge at the halfway point implying it is only $C^0$ continuos. For path 3, we notice that there is an sudden change in velocity around the halfway point showing that it is $C^1$ continuos.

#subpar.grid(
  figure(image("figs/E1.png")),
  figure(image("figs/E2.png")),
  grid.cell(figure(image("figs/E3.png", width: 50%)), colspan: 2),
  columns: (1fr, 1fr),
  caption: [Tracking Error for Each Path],
  label: <errors>,
)

Looking at the tracking error for each path, we see that both paths 1 and 2 have relatively comparable errors while path 3 almost has double the error of path 1 and it's maximum is reached at the halfway point. We also notice that the error for path 1 "spikes" suddenly, observing the robot at that point reveals that the arm is close to the robot body. This sudden spike might be attributed to the fact that the joint has to rotate at a wider angle to keep up with our desired path.

#subpar.grid(
  figure(image("figs/T1.png")),
  figure(image("figs/T2.png")),
  grid.cell(figure(image("figs/T3.png", width: 50%)), colspan: 2),
  columns: (1fr, 1fr),
  caption: [Force Quantities for Each Path],
  label: <forces>,
)

We note that we plot the curve of $f_3+19.6$ instead of $f_3$ as it remains relatively constant at that quantity as it constantly need to balance out the force exerted by the weight of the end-effector and the slight increases in height requested from it, which is why it sees very little fluctuations.

Observing the torque curve $tau_2$ for path 1 confirms the issue of being closer to the robot as the forces on that joint have to be higher to compensate for the wide change in angle. Perhaps this is why path 1 is not the optimal path despite being the shortest. We also notice a spike in torque in path 3 showing the sudden change of acceleration requested from the robot.

= Conclusion

In this project, we were able to analyze the performance of 3 different paths under the same control scheme. We concluded that path 2 (a linear path that remains away from the robot body) showed the best performance metrics and is most suitable under the same control scheme (Resolved Rate Control + PID), however our experiments are not encompassing enough to make a generalization for other robots, or even for different control schemes under the same configurations. 

This project helps law a solid foundation to begin the analysis of energy consumption across multiple configurations and control laws, though the methodology for a more general solution would have to be different and consider a wider range of paths and control parameters. For example, returning to #ref(<forces>) we could investigate higher degrees of smoothness to avoid sudden torque spikes. Perhaps a $C^2$ path might yield a better result?

#bibliography("sources.bib")