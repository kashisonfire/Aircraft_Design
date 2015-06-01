%% C.G. Calculation
% Run Main.m before running this code
%clc;
close all;
Length_VT = fusl+2;
%% Fuel Volume
K = tcavg+0.03;
xmax = b/2-fusd/2;
%150
MLE = 0.6538;
M10 = 0.6327;
M25 = 0.6009;
M70 = 0.5055;
MTE = 0.4420;
%220
% MLE = 40.97/69.42;
% M10 = 39.66/69.42;
% M25 = 37.68/69.42;
% M70 = 31.79/69.42;
% MTE = 27.86/69.42;
CR = 2*Sref/b/1.35;
CT = 0.35*CR;
mac = 2/3*(CR+CT-CR*CT/(CR+CT));
CRexp = 16.18;
Fx = @(x) K.*(MTE.*x+CRexp-MLE.*x).*(M70.*x+0.70*CRexp-M10.*x-0.10*CRexp);
FuelVolume = integral(Fx,0,xmax);

%% Center of Gravity for Fuel
syms y
Fz = K*(MTE*z+CRexp-MLE*z)*(M70*z+0.70*CRexp-M10*z-0.10*CRexp);
x_fuel = double(solve(int(Fz,0,y) == int(Fz,y,xmax), y));
y_fuel = .25*CR+x_fuel*tand(SA);

%% Weights
W_w = Ww*WTO^1.195;
W_fus = Wfus*WTO^0.235;
W_lg = 0.040*WTO;
W_np = 0.0555/WT*WTO;
W_ts = (0.17)*W_w;
W_fe = WFE+0.035*WTO;
W_pl = WPL;
W_pp = WTO/3.58/WT;
W_f = Wf2Wto*1.0275*WTO;


%% Center of Gravity (payload and fuel)
x_mac = b/6*(CR+2*CT)/(CR+CT);
y_mac = .25*CR+x_mac*tand(SA);
wing = 30;
for j = 1:50
    CGfus = 0.40*fusl*W_fus; %
    CGts = Length_VT*W_ts; %
    Dfan = 95.6/12*sqrt(Tengine/TSLT);
    Leng = 128.2/95.6*Dfan;
    Lnac = 1.1*(0.7*Dfan+Leng);
    CGpl = 0.52*fusl*W_pl; %
    CGw = W_w*(y_mac+wing); %
    CGfuel = W_f*(wing+y_fuel); %
    CGnp = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_np; %
    CGeng = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_pp; %
    CG = (CGfus+CGts+CGnp+CGeng+CGw+CGpl+CGfuel)/(WTO-W_fe-W_lg);
    COG_pf(j,:) = [wing CG];
    wing = wing+1;
end;

%% Center of Gravity (no payload)
wing = 50;
for j = 1:50
    CGfus = 0.40*fusl*W_fus; %
    CGts = Length_VT*W_ts; %
    Dfan = 95.6/12*sqrt(Tengine/TSLT);
    Leng = 128.2/95.6*Dfan;
    Lnac = 1.1*(0.7*Dfan+Leng);
    CGpl = 0; %
    CGw = W_w*(y_mac+wing); %
    CGfuel = W_f*(wing+y_fuel); %
    CGnp = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_np; %
    CGeng = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_pp; %
    CG = (CGfus+CGts+CGnp+CGeng+CGw+CGpl+CGfuel)/(WTO-W_pl-W_fe-W_lg);
    COG_np(j,:) = [wing CG];
    wing = wing+1;
end;


%% Center of Gravity (no fuel)
wing = 30;
for j = 1:50
    CGfus = 0.40*fusl*W_fus; %
    CGts = Length_VT*W_ts; %
    Dfan = 95.6/12*sqrt(Tengine/TSLT);
    Leng = 128.2/95.6*Dfan;
    Lnac = 1.1*(0.7*Dfan+Leng);
    CGpl = 0.52*fusl*W_pl; %
    CGw = W_w*(y_mac+wing); %
    CGfuel = 0; %
    CGnp = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_np; %
    CGeng = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_pp; %
    CG = (CGfus+CGts+CGnp+CGeng+CGw+CGpl+CGfuel)/(WTO-W_f-W_fe-W_lg);
    COG_nf(j,:) = [wing CG];
    wing = wing+1;
end;


%% Center of Gravity (no payload and fuel)
wing = 30;
for j = 1:50
    CGfus = 0.40*fusl*W_fus; %
    CGts = Length_VT*W_ts; %
    Dfan = 95.6/12*sqrt(Tengine/TSLT);
    Leng = 128.2/95.6*Dfan;
    Lnac = 1.1*(0.7*Dfan+Leng);
    CGpl = 0; %
    CGw = W_w*(y_mac+wing); %
    CGfuel = 0; %
    CGnp = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_np; %
    CGeng = (0.40*Lnac+wing+MLE*(.33*b/2+fusd/2)-2*1.1*Dfan)*W_pp; %
    CG = (CGfus+CGts+CGnp+CGeng+CGw+CGpl+CGfuel)/(WTO-W_pl-W_f-W_fe-W_lg);
    COG_npf(j,:) = [wing CG];
    wing = wing+1;
end;

mac25 = COG_npf(:,1)+y_mac;
mac15 = COG_npf(:,1)+y_mac-0.1*mac;
mac35 = COG_npf(:,1)+y_mac+0.1*mac;

figure(1)
plot(COG_pf(:,1),mac25,COG_pf(:,1),mac15,'--',COG_pf(:,1),mac35,'--',COG_pf(:,1),COG_pf(:,2),COG_nf(:,1),COG_nf(:,2),COG_np(:,1),COG_np(:,2),COG_npf(:,1),COG_npf(:,2))
legend mac25 mac15 mac35 both nofuel nopayload noboth
grid on
axis([40 60 55 80])
xlabel('Wing Position (ft)')
ylabel('C.G. (ft)')
% figure(2)
% plot(COG_pf(:,1),mac25,COG_pf(:,1),COG_pf(:,2),COG_nf(:,1),COG_nf(:,2),COG_np(:,1),COG_np(:,2),COG_npf(:,1),COG_npf(:,2))
% legend mac25 both nofuel nopayload noboth
% grid on
% axis([45 70 60 90])
