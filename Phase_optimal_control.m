clear;clc;close all;
%% Main file
% Phase-based control of periodic fluid flows
% Aditya Nair, Kunihiko Taira, Bingni Brunton, Steven Brunton
% September 2019
% Also cite: Monga, B., Wilson, D., Matchen, T., & Moehlis, J. (2019). 
% Phase reduction and phase-based optimal control for biological systems: 
% a tutorial. Biological cybernetics, 113(1-2), 11-46.

%% Load data
load ('Data/PRC_cylinder_blowing.mat');% Cylinder flow (momentum injection)
%load ('Data/PRC_cylinder_rotate.mat'); % Cylinder flow (rotary control)
%load ('Data/PRC_aoa09_blowing.mat');   % Airfoil flow (momentum injection)

%% Set parameters
global PS_x;global grad_PS_x;global Thet;global omega;
PS_x      = Phase_Sensitivity;               % Phase-sensitivity func. (PSF)
grad_PS_x = grad_Phase_Sensitivity;          % Gradient of PSF
Thet      = Theta;                           % Phase co-ordinate
T         = 6.1;                             % Time period (cylinder flow) (T = 1.141 for airfoil)
omega     = 2*pi/T;                          % Angular frequency
dt        = 0.005;                           % Time step
frac      = 0.8;                             % Fraction of time period
%(0.8: -2\pi/5 phase shift; 1.2: 2\pi/5 phase shift)
T_star    = frac*T;                          % Altered time period
eps       = 0.00001;
lambda    = 0.00001;
max_iterations = 10;

%% Solving BVP using Newton Iteration

for j = 1:max_iterations
    lambdap     = lambda + eps;yp0         = [0;lambdap];
    lambdam     = lambda - eps;ym0         = [0;lambdam];
    [~, yp]     = ode45(@euler_lagrange,[0 T_star],yp0);
    [~, ym]     = ode45(@euler_lagrange,[0 T_star],ym0);
    J           = (yp(end,1)-ym(end,1))/(2*eps)';
    y0          = [0;lambda];
    [~, y]      = ode45(@euler_lagrange,[0 T_star],y0);
    F           = (y(end,1)-2*pi);
    dF          = J;dlambda     = dF\(-F);
    lambda      = lambda + dlambda;
end
tspan           = 0:dt:T_star;
y0              = [0;lambda];
[t, y]          = ode45(@euler_lagrange,tspan,y0);
u               = y(:,2).*interp1(Thet,PS_x,wrapTo2Pi(y(:,1)))./(2);
theta           = y(:,1);
lambda          = y(:,2);

%% Euler-Lagrange equations

function df=euler_lagrange(~,y)
global PS_x;global grad_PS_x;global Thet;global omega;
df=zeros(2,1);
u    = y(2)*interp1(Thet,PS_x,wrapTo2Pi(y(1)))/(2);
df(1)= omega+interp1(Thet,PS_x,wrapTo2Pi(y(1)))*u;
df(2)=-u*y(2)*interp1(Thet,grad_PS_x,wrapTo2Pi(y(1)));
end
