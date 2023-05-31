close all 
clear all

%% OPERATING CONDITIONS
omega = 2*pi*6000/60; % rotor angular speed (rad/s)
rho = 1.2; % fluid density (kg/m3)

%% ROTOR GEOMETRY
Nb = 2; % number of blades
c = 0.025; % chord length (m)
r0 = 0.0125; % radial location of blade root (m)
R = 0.125; % rotor radius (m)
D = 2*R; % rotor diameter
sig = Nb*c/pi/R; % rotor solidity

%% AIRFOIL
Cla = 2*pi; % lift slope (/rad)

%% DISCRETIZE BLADE
N = 100; % number of blade elements
dr = (R-r0)/R/N; % non-dimensional element width
r = [r0/R:dr:1]+dr/2; % non-dimensional radial location

%% PITCH DISTRIBUTION
theta = ones(1,length(r)).*10*pi/180; % constant pitch

%% COMPUTE THRUST WITHOUT TIP LOSSES
lambda = (sig*Cla/16)*(sqrt(1+32.*theta.*r/sig/Cla)-1); 
dCt = (sig*Cla/2)*(theta.*r.^2-lambda.*r).*dr; % element thrust coefficient
Ct = sum(dCt); % total thrust coefficient
T = Ct*rho*(pi*R^2)*(omega*R)^2 % total thrust
figure(1); plot(r,dCt,'ko'); hold on;
phi = lambda./r; % inflow angle
figure(2); plot(r,(theta-phi).*180/pi,'ko'); hold on;

%% ADD TIP LOSSES
for i=1:10 % requires to be solved iteratively (typically converges after 4 or 5 iterations)
    phi = lambda./r;
    f = (Nb/2)*(1-r)./r./phi;
    F = (2/pi)*acos(exp(-f));
    lambda = (sig*Cla/16./F).*(sqrt(1+32.*F.*theta.*r/sig/Cla)-1);
    dCt = (sig*Cla/2)*(theta.*r.^2-lambda.*r).*dr;
    Ct = sum(dCt);
    T = Ct*rho*(pi*R^2)*(omega*R)^2
end
figure(1); plot(r,dCt,'ro');
figure(2); plot(r,(theta-phi).*180/pi,'ro');

%% FIGURES
figure(1); xlim([0 1]);
set(gca,'fontsize',18);
x=xlabel('r/R','fontSize', 24,'Interpreter','latex');
ylabel('$\Delta C_t$','fontSize', 24,'Interpreter','latex');
legend('no tip loss','with tip losses'); legend('location','northwest');
figure(2); xlim([0 1]);
set(gca,'fontsize',18);
x=xlabel('r/R','fontSize', 24,'Interpreter','latex');
ylabel('$\alpha$ ($^{\circ}$)','fontSize', 24,'Interpreter','latex');
legend('no tip loss','with tip losses'); legend('location','northwest');

