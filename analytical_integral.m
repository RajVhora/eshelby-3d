%Analytical solution of Eshelby problem for Homogeneous inclusion using
%Elliptical Integral.
clc
close all
clear all
a1 = 1e20;
a2 = 50.000000001;
a3 = 50;
theta = asin(sqrt(1-a3^2/a1^2));
k = ((a1^2-a2^2)/(a1^2-a3^2));
F = ellipticF(theta,k);
E = ellipticE(theta,k);
I1 = 4*pi*a1*a2*a3*(F-E)/(a1^2-a2^2)/sqrt(a1^2-a3^2);
I2 = 4*pi*a1*a2*a3*(a2*sqrt(a1^2-a3^2)/a1/a3-E)/(a2^2-a3^2)/sqrt(a1^2-a3^2);
I3 = 4*pi - I1-I2;
I12 = (I2-I1)/(a1^2-a2^2);
syms I11 I13
eqn1 = 3*I11 +I13 == 4*pi/a1^2 -I12;
eqn2 = 3*a1^2*I11 +a3^2*I13 == 3*I1 -a2^2*I12;
[A,B] = equationsToMatrix([eqn1,eqn2],[I11,I13]);
 p = linsolve(A,B);
 I11 = p(1);
 I13 = p(2);

%%%%%%%%%%%%%%%%%%%%%%
I23 = (I3-I2)/(a2^2-a3^2);
syms I22 I21
eqn1 = 3*I22 +I21 == 4*pi/a2^2 -I23;
eqn2 = 3*a2^2*I22 +a1^2*I21 == 3*I2 -a3^2*I23;
[A,B] = equationsToMatrix([eqn1,eqn2],[I22,I21]);
q= linsolve(A,B);
I22 = q(1);
I21 = q(2);

%%%%%%%%%%%%%%%%%%%
I31 = (I1-I3)/(a3^2-a1^2);
syms I33 I32
eqn1 = 3*I33 +I32 == 4*pi/a3^2 -I31;
eqn2 = 3*a3^2*I33 +a2^2*I32 == 3*I3 -a1^2*I31;
[A,B] = equationsToMatrix([eqn1,eqn2],[I33,I32]);
r = linsolve(A,B);
I33 = r(1);
I32 = r(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 200;
nu = 3/10;
epsilon11 = 0.00;
epsilon22 = 0.01;
epsilon33 = 0.01;
epsilon23 = 0.00;

x = (a2^2/8/pi/(1-nu)*((1-nu)/(1-2*nu)*3*I22 + nu/(1-2*nu)*(I32 + I12))....
    + (1-2*nu)/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I2 -nu/(1-2*nu)*(I3+I1))....
    - (1-nu)/(1-2*nu))*epsilon22;

y = (a3^2/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I23 + nu/(1-2*nu)*(3*I33 + I13))....
    - (1-2*nu)/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I2 -nu/(1-2*nu)*(I3-I1))....
    - nu/(1-2*nu))*epsilon33;
z = (a1^2/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I21 + nu/(1-2*nu)*(3*I11 + I31))....
    - (1-2*nu)/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I2 - nu/(1-2*nu)*(I1-I3))....
    - nu/(1-2*nu))*epsilon11;


sigma22 = x + y + z;

sigma23 = ((a2^2+a3^3)/8/pi/(1-nu)*I23 + (1-2*nu)/8/pi/(1-nu)*(I2+I3) -1)*epsilon23;

x1 = (a3^2/8/pi/(1-nu)*((1-nu)/(1-2*nu)*3*I33 + nu/(1-2*nu)*(I13 + I23))....
    + (1-2*nu)/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I3 -nu/(1-2*nu)*(I1+I2))....
    - (1-nu)/(1-2*nu))*epsilon33;

y1 = (a1^2/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I31 + nu/(1-2*nu)*(3*I11 + I21))....
    - (1-2*nu)/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I3 -nu/(1-2*nu)*(I1-I2))....
    - nu/(1-2*nu))*epsilon11;
z1 = (a2^2/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I32 + nu/(1-2*nu)*(3*I22 + I12))....
    - (1-2*nu)/8/pi/(1-nu)*((1-nu)/(1-2*nu)*I3 - nu/(1-2*nu)*(I2-I1))....
    - nu/(1-2*nu))*epsilon22;

sigma33 = x1 + y1 + z1;

Sigmaxx = double(sigma22)*2/0.01
Sigmayy = double(sigma33)*2/0.01
Sigmaxy = double(sigma23)*2/0.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
