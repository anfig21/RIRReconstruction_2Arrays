%% KERNEL FUNCTIONS
clear, clc, close all

%% Initial parameters
c = 343;
f = 1e3;
k = 2*pi*f/c;
Delta_r = linspace(0,10,1000);

%% Kernel
K_Sph = sinc(k*Delta_r/pi);
K_Cyl = sin(k*Delta_r)./sqrt(k*Delta_r);

%% Plot
figure, hold on
plot(Delta_r,K_Sph), grid on
plot(Delta_r,K_Cyl)

legend('Spherical','Cylindrical')