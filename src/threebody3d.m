%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function threebody3d
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
This function solves the inertial three-body problem in three dimensions
numerically using the RKF4(5) method.

G                   -universal gravitational constant (km^3/kg s^2)
mE, mS1, mS2        -the masses of the two bodies (Earth and satellite)(kg)
m1                  -total mass (mE + mS1)
m2                  -total mass (mE + mS2)
t0                  -initial time (s)
tf                  -final time (s)
RS1_0, VS1_0        -3 by 1 column vector containing the components of the
                    initial position (km) and velocity (km/s) of mS1
RS2_0, VS2_0        -3 by 1 column vector containing the componenets of the
                    initial position (km) and velocity (km/s) of mS2
y0                  -12 by 1 column vector containing the initial values of
                    the state vectors of the two bodies:
                    [RS1_0; RS2_0; VS1_0; VS2_0]
t                   -column vector of the times at which the solution is
                    found
XS1, YS1, ZS1       -column vectors containing the X,Y,Z coordinates (km)
                    of mS1 at the times in t
XS2, YS2, ZS2       -column vectors containing the X,Y,Z coordinates (km)
                    of mS2 at the times in t
VXS1, VYS1, VZS1    -column vectors containing Vx, Vy, Vz coordinates
                    (km/s) of mS1 at the times in t
VXS2, VYS2, VZS2    -column vectors containing Vx, Vy, Vz coordinates
                    (km/s) of mS2 at the times in t
y                   -a matrix whose 12 columns are, respectively,
                    XS1, YS1, ZS1, XS2, YS2, ZS2, VXS1, VYS1, VZS1, VXS2,
                    VYS2, VZS2

User M-functions required: rkf45
User subfunctions required: rates, output
%}
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc; clear all; close all

G                   = 6.67259e-20;

%...Input data:
mE                  = 100.e26;
mS1                 = 20.e26;
mS2                 = 1.e26;
t0                  = 0;
tf                  = 100000;
timespan            = [t0 tf];

RS1_0               = [500000;    0;  0];
RS2_0               = [600000;    0;  0];
VS1_0               = [0;     45;  5];
VS2_0               = [0;     72;  10];
%...End input data

y0                  = [RS1_0; RS2_0; VS1_0; VS2_0];

%...Integrate the equations of motion:
[t, y]              = rkf45(@rates, timespan, y0);

%...Output the results:
output

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dydt = rates(t,y)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
This function calculates the accelerations for mS1 and mS2

t                   -time (s)
y                   -column vector containing the postion and velocity
                    vectors of the system at time t
RS1, RS2            -position vectors of mS1 and mS2
VS1, VS2            -velocity vectors of mS1 and mS2
r1                  -magnitude of the relative position from mS1 to mE
r2                  -magnitude of the relative position from mS2 to mE
RS, rS              -vector and magnitude of the relative position from mS1 to mS2
AS1, AS2            -acceleration vectors of mS1 and mS2
dydt                -column vectors containing the velocity and
                    acceleration vectors of the sytem at time t
%}
%----------------------------------

RS1                 = [y(1); y(2); y(3)];
RS2                 = [y(4); y(5); y(6)];
VS1                 = [y(7); y(8); y(9)];
VS2                 = [y(10); y(11); y(12)];

m1                  = (mS1 + mE);
m2                  = (mS2 + mE);
r1                  = norm(RS1);
r2                  = norm(RS2);
RS                  = (RS2 - RS1);
rs                  = norm(RS);

AS1                 = -(G*m1*RS1/r1^3) + (G*mS2*RS/rs^3);
AS2                 = -(G*m2*RS2/r2^3) - (G*mS1*RS/rs^3);

dydt = [VS1; VS2; AS1; AS2];

end %rates
%-------------------------------------

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function output
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
This function plots the motion of mS1 and mS2 relative to mE at the origin
%}

%...Extract the particle trajectories:
XS1 = y(:,1);
YS1 = y(:,2);
ZS1 = y(:,3);
XS2 = y(:,4);
YS2 = y(:,5);
ZS2 = y(:,6);

%...Plot the trajectories

figure (1)
title('Motion of 2 satellites relative to larger body')
hold on
plot3(XS1, YS1, ZS1, '-r')
plot3(XS2, YS2, ZS2, '-b')
text(0, 0, 0, 'o')
%axis('equal')
view([2,45])
grid on
%axis equal
xlim([-1500000 1500000])
ylim([-1500000 1500000])
zlim([-1500000 1500000])
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')
end %output

end %threebody3d
