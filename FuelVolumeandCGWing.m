% Run Main.m before running this code
clc;
close all;

%% Fuel Volume
K = tcavg+0.03;
xmax = b/2-fusd/2;
MLE = 0.6538;
M10 = 0.6327;
M25 = 0.6009;
M70 = 0.5055;
MTE = 0.4420;
CR = 2*Sref/b/1.35;
CT = 0.35*CR;
mac = 2/3*(CR+CT-CR*CT/(CR+CT));
CRexp = 16.38;
Fx = @(x) K.*(MTE.*x+CRexp-MLE.*x).*(M70.*x+0.70*CRexp-M10.*x-0.10*CRexp);
FuelVolume = integral(Fx,0,xmax);

%% Center of Gravity for Fuel
length = 2000;
for i = 1:length
    V(i) = integral(Fx,xmax/length*(i-1),xmax/length*i);
    pos(i) = xmax/length*(i-1)+.5*xmax/length;
    Vx(i) = pos(i)*V(i);
end;
xfuel = sum(Vx)/FuelVolume*M25;
