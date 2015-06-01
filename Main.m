%% MAE 159 Aircraft Performance
% This program will allow the user to test a number of different aircraft
% designs by changing design specifications, wing and engine parameters,
% and the internal configuration of the airplane.
addpath(genpath('C:\Users\Akash\Google Drive\MAE 159'))
tic;
clear all;
close all;
format long;

%% Iterations, Varying Factors, and Output
VaryAR = 'no'; % Vary Aspect Ratio? (only vary one of these)
VarySA = 'no'; % Vary Sweep Back?
VaryFC= 'no'; %varyfuelcost
InitialSA = 31; % Initial Value
InitialAR = 9.1; % Initial Value
InitialFC = 0.28; % Initial Value
Increment = .5; % Increment between values
FinalValue = 35; % Final Value of SA or AR
output = 0; % 1 = Output to Excel, 2 = Nothing
filename = '150PAXsup.xlsx'; % File name
sheet = 'compareconsup'; % Excel Sheet Desired
Row = 55; % Row location output is desired
addTitle = 0; % If you want the output to include a Title

%% Design Specifications
NoPAX = 150;      % Number of Passengers
Wcargo = 3000;    % Weight of Cargo in pounds (10 pounds/ft^3)
Range = 3500;     % Nautical miles (Still air)
TOFL = 7500;      % Takeoff field length in feet (sea-level, hot day(84 F)
Vapp = 135;       % Landing approach speed
MaxFuel = 0.25;   % Use 0.75 for Liebeck, and 0.25 for our Design Spec
class3 = 1.0;     % Use 1.1 for Liebeck, and 1.0 for our Design Spec
MCruise = 0.80;   % Cruise Mach number
Altitude = 35000; % Initial Cruise Altitude
FuelCost = InitialFC;  % 28 cents for Domestic, 41 cents for International in 1975

%% Wing Parameters
wtype = 'sup';         % Type of Wing (enter con or sup)
SA = InitialSA;        % Sweepback Angle
AR = InitialAR;        % Aspect Ratio
tratio = 0.35;         % Taper Ratio

%% Engine Parameters
etype = 'JT9D';        % Engine Type (JT8D or JT9D)
NoEng = 2;             % Number of engines (2,3,4)
Location = 0;          % Location of Engine (0 for Wing, 1 for Fuselage)

%% Interior Configuration Parameters
NoAisles = 1;          % Number of aisles (1 or 2)
NoAbreast = 6;         % Number of seats abreast (4-8)

%% Advance Wing Technology
advtech = 1;           % Advanced Technology of Materials
comp = 1;              % Composite Material
alloy = 0;             % Aluminum/Lithium Material

%% Advance Engine Technology
AdvancedEngineTec = 1; % Advanced Technology of Engine

%% Tolerance
tolCL = 1e-7;          % Accuracy of CL to CL_IC
tolWTO = 1e-3;         % Accuracy of Take Off Weight
tolRange = 20;         % Accuracy of Range to All Out Range
tolWT = 1e-3;          % Accuracy of Weight2Thrust Increase with Gradient
rateofinc = 1e-5;      % Rate at which the Wf2Wto will increase (1e-5 is best)

if strcmp(VaryAR,'yes')
    length2 = (FinalValue-InitialAR)/Increment+1;
elseif strcmp(VarySA,'yes')
    length2 = (FinalValue-InitialSA)/Increment+1;
elseif strcmp(VaryFC,'yes')
    length2 = (FinalValue-InitialFC)/Increment+1;
else
    length2 = 1;
end;
for q=1:length2
    %% Initiate Total Range/Range All Out Loop
    fail = 0;
    checkc = 0;
    Vcruise = MCruise*(0.592484*sqrt(1.4*1716.49*(394.4)));
    Rao = Range+200+.75*Vcruise;
    fig4 = csvread('fig4Data.txt');
    Wf2WtoJT8D = interp1(fig4(:,1),fig4(:,2),Rao);
    Wf2Wto = zeros(2,1);
    if strcmp(etype,'JT8D') == 1
        Wf2Wto(1) = Wf2WtoJT8D;
    else
        Wf2Wto(1) = Wf2WtoJT8D*0.61/0.78;
    end;
    maxit = 30;
    Wf2WtoConverge = 0;
    GradConverge = 0;
    iteration = 1;
    storage = 9999999;
    for z=1:maxit
        %% Wing Design Calculation
        start = 0;
        if Wf2WtoConverge ~= 1 && fail~=1
            CL = zeros(1,2);
            CL(1) = 0.5;
            itmax = 1000;
            for j = 1:itmax
                if strcmp(wtype,'con') == 1
                    DeltaMDIV = -0.2412*(CL(j))^2-0.0671*CL(j)+0.1099;
                elseif strcmp(wtype,'sup') == 1
                    DeltaMDIV = 508.6*(CL(j))^6 - 1685.1*(CL(j))^5 + 2309.7*(CL(j))^4 - 1675.3*(CL(j))^3 + 677.37*(CL(j))^2 - 144.61*(CL(j)) + 12.741;
                end;
                MDIV = MCruise + 0.004 - DeltaMDIV;
                if strcmp(wtype,'con') == 1
                    if SA >= 0 && SA < 10
                        tc1 = -0.6363*MDIV + 0.5742;
                        tc2 = -0.6171*MDIV + 0.5643;
                        tcavg = tc1 - (tc1-tc2)/(0-10)*(0-SA);
                    elseif SA >= 10 && SA < 15
                        tc1 = -0.6171*MDIV + 0.5643;
                        tc2 = -0.5939*MDIV + 0.5521;
                        tcavg = tc1 - (tc1-tc2)/(10-15)*(10-SA);
                    elseif SA >= 15 && SA < 20
                        tc1 = -0.5939*MDIV + 0.5521;
                        tc2 = -0.5662*MDIV + 0.5373;
                        tcavg = tc1 - (tc1-tc2)/(15-20)*(15-SA);
                    elseif SA >= 20 && SA < 25
                        tc1 = -0.5662*MDIV + 0.5373;
                        tc2 = -0.5331*MDIV + 0.5196;
                        tcavg = tc1 - (tc1-tc2)/(20-25)*(20-SA);
                    elseif SA >= 25 && SA < 30
                        tc1 = -0.5331*MDIV + 0.5196;
                        tc2 = -0.5055*MDIV + 0.5062;
                        tcavg = tc1 - (tc1-tc2)/(25-30)*(25-SA);
                    elseif SA >= 30 && SA < 35
                        tc1 = -0.5055*MDIV + 0.5062;
                        tc2 = -0.4699*MDIV + 0.4873;
                        tcavg = tc1 - (tc1-tc2)/(30-35)*(30-SA);
                    elseif SA >= 35 && SA <= 40
                        tc1 = -0.4699*MDIV + 0.4873;
                        tc2 = -0.4301*MDIV + 0.4654;
                        tcavg = tc1 - (tc1-tc2)/(35-40)*(35-SA);
                    end;
                elseif strcmp(wtype,'sup') == 1
                    if (0<=SA) && (SA<5)
                        tc1 = 0.0208*MDIV^-5.653;
                        tc2 = 0.0216*MDIV^-5.61;
                        tcavg = tc1+((SA/5)*(tc2-tc1));
                    end;
                    if (5<=SA) && (SA<10)
                        tc1 = 0.0216*MDIV^-5.61;
                        tc2 = 0.0242*MDIV^-5.406;
                        tcavg = tc1+((SA/5)*(tc2-tc1));
                    end;
                    if (10<=SA) && (SA<15)
                        tc1 = 0.0242*MDIV^-5.406;
                        tc2 = 0.0282*MDIV^-5.177;
                        tcavg= tc1+((SA-10)/5)*(tc2-tc1);
                    end;
                    if (15<=SA) && (SA<20)
                        tc1 = 0.0282*MDIV^-5.177;
                        tc2 = 0.0329*MDIV^-5.006;
                        tcavg = tc1+((SA-15)/5)*(tc2-tc1);
                    end;
                    if (20<=SA) && (SA<25)
                        tc1 = 0.0329*MDIV^-5.006;
                        tc2 = 0.0368*MDIV^-5.107;
                        tcavg= tc1+((SA-20)/5)*(tc2-tc1);
                    end;
                    if (25<=SA) && (SA<30)
                        tc1 = 0.0368*MDIV^-5.107;
                        tc2 = 0.0394*MDIV^-5.577;
                        tcavg= tc1+((SA-25)/5)*(tc2-tc1);
                    end;
                    if (30<=SA) && (SA<35)
                        tc1 = 0.0394*MDIV^-5.577;
                        tc2 = 0.0391*MDIV^-6.854;
                        tcavg= tc1+((SA-30)/5)*(tc2-tc1);
                    end;
                    if (35<=SA) && (SA<=40)
                        tc1 = 0.0391*MDIV^-6.854;
                        tc2 = 0.0320*MDIV^-10.75;
                        tcavg= tc1+((SA-35)/5)*(tc2-tc1);
                    end;
                end;
                Relation1 = (cosd(SA))^2*(tcavg)^2*AR;
                CLmaxto = -7489*Relation1^5 + 6045.9*Relation1^4 - 1742.2*Relation1^3 + 189.44*Relation1^2 + 1.0242*Relation1 + 1.3951;
                CLmaxldg = -8066.4*Relation1^5 + 4998.9*Relation1^4 - 1029.7*Relation1^3 + 49.013*Relation1^2 + 11.275*Relation1 + 2.0813;
                WSldg = (Vapp/1.3)^2*.953.*CLmaxldg./296;
                WSto = WSldg./(1-MaxFuel*Wf2Wto(z));
                WSIC = .965.*WSto;
                CL(j+1) = WSIC./(1481*.2360*.82^2);
                diff = abs(CL(j+1)-CL(j));
                if diff < tolCL
                    break;
                end;
            end;
            sizeCL = size(CL);
            CL = CL(sizeCL(1,2));
            %% Takeoff Field Length
            if NoEng == 2
                WT07VLO = (-2.881e-07*TOFL^2 + 3.262e-02*TOFL - 2.329e+01)*(.953*CLmaxto/WSto);
            end;
            if NoEng == 3
                WT07VLO = (-3.0179e-07*TOFL^2 + 3.5921e-02*TOFL - 2.1486e+01)*(.953*CLmaxto/WSto);
            end;
            if NoEng == 4
                WT07VLO = (-3.5052e-07*TOFL^2 + 3.7488e-02*TOFL - 1.3539e+01)*(.953*CLmaxto/WSto);
            end;
            VLO = 1.2*(296*WSto/(0.953*CLmaxto))^.5;
            MLO = VLO/(0.592484*sqrt(1.4*1716.49*(84+459.67)));
            M07LO = 0.7*MLO;
            if strcmp(etype,'JT8D') == 1
                TSLT = 14500;
                Tm = 7606.2*M07LO^2-9185.4*M07LO+14401;
                WT = WT07VLO*Tm/TSLT;
            elseif strcmp(etype,'JT9D') == 1
                TSLT = 45500;
                figSLJT9D = csvread('SLDTOT.txt');
                Tm = 37470*M07LO^2 - 47380*M07LO + 45368;
                WT = WT07VLO*Tm/TSLT;
            end;
        end;
        ClimbCheck = 0;
        while start==0
            if fail == 1
                break;
            end;
            %% Weight Estimation
            if ClimbCheck == 0
                syms WTO
                Ww = .00945*AR^.8*(1+tratio)^0.25*1.01*3.75^.5/((tcavg+0.03)^.4*cosd(SA)*(WSto)^0.695);
                fusl = (3.76*NoPAX/NoAbreast+33.2)*class3;
                fusd = (1.75*NoAbreast+1.58*NoAisles+1)*class3;
                Wfus = .6727*11.5*fusl^.6*fusd^.72*3.75^.3;
                Wlg = 0.04;
                Wnp = 0.0555/WT;
                WwTS = (0.17+0.08/NoEng*Location)*Ww+Ww;
                WPP = 1/(WT*3.58);
                Wfuel = Wf2Wto(z)*1.0275;
                WPL = 215*NoPAX+Wcargo;
                WFE = 132*NoPAX+300*NoEng+260*2+170*(ceil(NoPAX/50));
                itmax = 100;
                WTO = zeros(1,2);
                WTO(1) = 200000;
                WTOal = zeros(1,2);
                WTOal(1) = 200000;
                if AdvancedEngineTec == 1
                    WPP = 1.1*WPP;
                end;
                for n = 1:itmax
                    WTOal(n+1) = WwTS*(WTOal(n))^1.195+Wfus*(WTOal(n))^0.235+(Wnp+Wlg+WPP+Wfuel+0.035)*(WTOal(n))+WFE+WPL;
                    diff = abs(WTOal(n+1)-WTOal(n));
                    if diff < tolWTO
                        break;
                    end;
                end;
                WTOal = WTOal(length(WTO));
                if advtech == 1 && comp == 1
                    WwTS = WwTS*0.7;
                    Wfus = Wfus*0.85;
                    WFE = WFE*0.9;
                    Wnp = 0.8*Wnp;
                end;
                if advtech == 1 && alloy == 1
                    WwTS = WwTS*0.94;
                    Wfus = Wfus*0.94;
                end;
                for n = 1:itmax
                    WTO(n+1) = WwTS*(WTO(n))^1.195+Wfus*(WTO(n))^0.235+(Wnp+Wlg+WPP+Wfuel+0.035)*(WTO(n))+WFE+WPL;
                    diff = abs(WTO(n+1)-WTO(n));
                    if diff < tolWTO
                        break;
                    end;
                end;
                WTO = WTO(length(WTO));
                
                %% Important Aircraft Parameters
                Sref = WTO/WSto;         % Wing Planform Area (ft^2)
                b = (AR*Sref)^.5;        % Wing Span (ft)
                %Sref = 1912.22;
                %b = 139.66;
                mac = Sref/b;            % Mean Aerodynamic Chord (ft)
                Thrust = WTO/WT;         % Thrust (lb)
                Tengine = Thrust/NoEng;  % Thrust per engine (lb)
                
                %% Drag
                % Skin friction approximated using Figure 11.2 From Shevell
                % Cf is a function of Reynolds Number
                % Density, Velocity, and Dynamic Viscosity obtained at
                % M_Cruise = 0.50 and Altitude = 30,000 feet
                
                RNperFT = 8.91e-4*.5*sqrt(1.4*1718*(459.67-47.83))/3.107e-7;
                fig112 = csvread('SkinFriction.txt');
                
                % Wing/Tail
                RNw = RNperFT*mac;
                natRNw = log(RNw);
                Cfw = 1e-3*exp(6.42329e-06*natRNw^6 - 6.12758e-04*natRNw^5 + 2.41827e-02*natRNw^4 - 5.05612e-01*natRNw^3 + 5.91688e+00*natRNw^2 - 3.70536e+01*natRNw + 1.00265e+02);
                SwetW = 2*(Sref-(0.85*2*mac/(1+tratio))*fusd)*1.02;
                if SA >= 0 && SA < 10
                    Kw = 267.45*tcavg^4 - 70.235*tcavg^3 + 9.7108*tcavg^2 + 1.5765*tcavg + 1.0008;
                elseif SA >= 10 && SA < 15
                    Kw1 = 267.45*tcavg^4 - 70.235*tcavg^3 + 9.7108*tcavg^2 + 1.5765*tcavg + 1.0008;
                    Kw2 = 28.578*tcavg^3 - 1.6311*tcavg^2 + 1.8873*tcavg + 0.9989;
                    Kw = Kw1 - (Kw1-Kw2)/(10-15)*(10-SA);
                elseif SA >= 15 && SA < 20
                    Kw1 = 28.578*tcavg^3 - 1.6311*tcavg^2 + 1.8873*tcavg + 0.9989;
                    Kw2 = 153.65*tcavg^4 - 30.71*tcavg^3 + 5.8145*tcavg^2 + 1.4239*tcavg + 1.003;
                    Kw = Kw1 - (Kw1-Kw2)/(15-20)*(15-SA);
                elseif SA >= 20 && SA < 25
                    Kw1 = 153.65*tcavg^4 - 30.71*tcavg^3 + 5.8145*tcavg^2 + 1.4239*tcavg + 1.003;
                    Kw2 = 282.52*tcavg^4 - 70.398*tcavg^3 + 9.0241*tcavg^2 + 1.3817*tcavg + 0.9999;
                    Kw = Kw1 - (Kw1-Kw2)/(20-25)*(20-SA);
                elseif SA >= 25 && SA < 30
                    Kw1 = 282.52*tcavg^4 - 70.398*tcavg^3 + 9.0241*tcavg^2 + 1.3817*tcavg + 0.9999;
                    Kw2 = 162.49*tcavg^4 - 29.447*tcavg^3 + 4.508*tcavg^2 + 1.4476*tcavg + 1.0008;
                    Kw = Kw1 - (Kw1-Kw2)/(25-30)*(25-SA);
                elseif SA >= 30 && SA < 35
                    Kw1 = 162.49*tcavg^4 - 29.447*tcavg^3 + 4.508*tcavg^2 + 1.4476*tcavg + 1.0008;
                    Kw2 = 212.41*tcavg^4 - 48.243*tcavg^3 + 6.7731*tcavg^2 + 1.2485*tcavg + 1.0015;
                    Kw = Kw1 - (Kw1-Kw2)/(30-35)*(30-SA);
                elseif SA >= 35 && SA <= 40
                    Kw1 = 212.41*tcavg^4 - 48.243*tcavg^3 + 6.7731*tcavg^2 + 1.2485*tcavg + 1.0015;
                    Kw2 = 194.6*tcavg^4 - 38.877*tcavg^3 + 5.1437*tcavg^2 + 1.2233*tcavg + 1.0009;
                    Kw = Kw1 - (Kw1-Kw2)/(35-40)*(35-SA);
                end
                fwing = Kw*Cfw*SwetW;
                ftail = (0.35+ 0.10/NoEng*Location)*fwing;
                
                % Fuselage
                RNf = RNperFT*fusl;
                natRNf = log(RNf);
                Cff = 1e-3*exp(6.42329e-06*natRNf^6 - 6.12758e-04*natRNf^5 + 2.41827e-02*natRNf^4 - 5.05612e-01*natRNf^3 + 5.91688e+00*natRNf^2 - 3.70536e+01*natRNf + 1.00265e+02);
                Swetf = 0.9*pi*fusl*fusd;
                LDratio = fusl/fusd;
                Kf = -0.00006644*LDratio^5 + 0.00269185*LDratio^4 - 0.04357548*LDratio^3 + 0.35480582*LDratio^2 - 1.49552252*LDratio + 3.86618516;
                ffus = Kf*Cff*Swetf;
                
                % Nacelles/Pylon
                SwetN = 2.1*(Tengine)^.5*NoEng;
                fnac = 1.25*Cfw*SwetN;
                fpylon = 0.20*fnac;
                
                ftotal = (fwing+ftail+ffus+fnac+fpylon)*1.06;
                CDo = ftotal/Sref;
                e = 1/(1.035+0.38*CDo*pi*AR);
                
                %% Climb
                Vcl = 1.3*12.9/(ftotal*e)^.25*((WTO-.0175*WTO)/(0.5702*b))^.5;
                Mcl = Vcl/(0.592484*sqrt(1.4*1718*(-12.26+459.67)));
                TR = 0.5702*ftotal*Vcl^2/296+94.1/(0.5702*e*Vcl^2)*(WTO/b)^2;
                if strcmp(etype,'JT9D') == 1
                    TA1 = -23053*Mcl^5 + 73875*Mcl^4 - 91469*Mcl^3 + 60201*Mcl^2 - 30922*Mcl + 27143;
                    c1 = 7.604e-02*Mcl^2 + 3.471e-01*Mcl + 3.354e-01;
                    TA2 = 157192*Mcl^6 - 583328*Mcl^5 + 856431*Mcl^4 - 634808*Mcl^3 + 252828*Mcl^2 - 54627*Mcl + 20289;
                    c2 = 0.3698*Mcl + 0.3417;
                    cClimb = (c1+c2)/2;
                    TAcl = (TA1+TA2)/2;
                end;
                if strcmp(etype,'JT8D') == 1
                    TAcl = 36179*Mcl^6+122034*Mcl^5-150591*Mcl^4+81423*Mcl^3-14946*Mcl^2-2643.2*Mcl+7472;
                    cClimb = -0.9887*Mcl^6+2.1185*Mcl^5-0.7239*Mcl^4-1.216*Mcl^3+0.9944*Mcl^2+0.0909*Mcl+0.5905;
                end;
                if AdvancedEngineTec == 1
                    cClimb = cClimb*0.9;
                end;
                TA = Tengine/TSLT*TAcl;
                RClimb = 101*(NoEng*TA-TR)*Vcl/WTO;
                TimeCl = Altitude/RClimb;
                RangeCl = Vcl*TimeCl/60;
                WfCl = NoEng*TA*cClimb*TimeCl/60;
                
                %% Range
                W_0 = WTO-WfCl;
                W_1 = (1-Wf2Wto(z))*WTO;
                Cl_cruise = (W_0+W_1)/(2*Sref)/(1481*0.2360*MCruise^2);
                Cd = Cl_cruise^2/(pi*AR*e)+CDo+0.001;
                L2D = Cl_cruise/Cd;
                Trr = (W_0+W_1)/2/L2D;
                TrJT9D = Trr*TSLT/Tengine/NoEng;
                if strcmp(etype,'JT8D') == 1
                    c = 9.553e-15*TrJT9D^5-1.814e-10*TrJT9D^4+1.376e-6*TrJT9D^3-5.218e-3*TrJT9D^2+9.884*TrJT9D-7484;
                end;
                if strcmp(etype,'JT9D') == 1
                    c = 3.3144e-17*TrJT9D^4 - 1.3542e-12*TrJT9D^3 + 2.2789e-08*TrJT9D^2 - 1.6998e-04*TrJT9D + 1.0587; % Fit equation from 35000 ft chart
                end;
                if AdvancedEngineTec == 1
                    c = c*0.9;
                end;
                R_cruise = Vcruise/c*L2D*log(W_0/W_1);
                RangeTotal = RangeCl+R_cruise;
                if z > 2
                    storage = diff2;
                end;
                diff2 = abs(RangeTotal-Rao);
                if storage < diff2 && checkc == 0
                    rateofinc = -1*rateofinc;
                    checkc = 1;
                end;
                Wf2Wto(z+1) = Wf2Wto(z)+rateofinc*diff2;
            elseif ClimbCheck == 1
                syms WTO
                Ww = .00945*AR^.8*(1+tratio)^0.25*1.01*3.75^.5/((tcavg+0.03)^.4*cosd(SA)*(WSto)^0.695);
                fusl = (3.76*NoPAX/NoAbreast+33.2)*class3;
                fusd = (1.75*NoAbreast+1.58*NoAisles+1)*class3;
                Wfus = .6727*11.5*fusl^.6*fusd^.72*3.75^.3;
                Wlg = 0.04;
                Wnp = 0.0555/WT;
                WwTS = (0.17+0.08/NoEng*Location)*Ww+Ww;
                WPP = 1/(WT*3.58);
                Wfuel = Wf2Wto(z)*1.0275;
                WPL = 215*NoPAX+Wcargo;
                WFE = 132*NoPAX+300*NoEng+260*2+170*(ceil(NoPAX/50));
                itmax = 100;
                WTO = zeros(1,2);
                WTO(1) = 200000;
                WTOal = zeros(1,2);
                WTOal(1) = 200000;
                if AdvancedEngineTec == 1
                    WPP = 1.1*WPP;
                end;
                for n = 1:itmax
                    WTOal(n+1) = WwTS*(WTOal(n))^1.195+Wfus*(WTOal(n))^0.235+(Wnp+Wlg+WPP+Wfuel+0.035)*(WTOal(n))+WFE+WPL;
                    diff = abs(WTOal(n+1)-WTOal(n));
                    if diff < tolWTO
                        break;
                    end;
                end;
                WTOal = WTOal(length(WTO));
                if advtech == 1 && comp == 1
                    WwTS = WwTS*0.7;
                    Wfus = Wfus*0.85;
                    WFE = WFE*0.9;
                    Wnp = 0.8*Wnp;
                end;
                if advtech == 1 && alloy == 1
                    WwTS = WwTS*0.94;
                    Wfus = Wfus*0.94;
                end;
                for n = 1:itmax
                    WTO(n+1) = WwTS*(WTO(n))^1.195+Wfus*(WTO(n))^0.235+(Wnp+Wlg+WPP+Wfuel+0.035)*(WTO(n))+WFE+WPL;
                    diff = abs(WTO(n+1)-WTO(n));
                    if diff < tolWTO
                        break;
                    end;
                end;
                WTO = WTO(length(WTO));
                
                %% Important Aircraft Parameters
                Sref = WTO/WSto;         % Wing Planform Area (ft^2)
                b = (AR*Sref)^.5;        % Wing Span (ft)
                %Sref = 1912.22;
                %b = 139.66;
                mac = Sref/b;            % Mean Aerodynamic Chord (ft)
                Thrust = WTO/WT;         % Thrust (lb)
                Tengine = Thrust/NoEng;  % Thrust per engine (lb)
                
                %% Drag
                % Skin friction approximated using Figure 11.2 From Shevell
                % Cf is a function of Reynolds Number
                % Density, Velocity, and Dynamic Viscosity obtained at
                % M_Cruise = 0.50 and Altitude = 30,000 feet
                
                RNperFT = 1.426*10^6;
                fig112 = csvread('SkinFriction.txt');
                
                % Wing/Tail
                RNw = RNperFT*mac;
                natRNw = log(RNw);
                Cfw = 1e-3*exp(6.42329e-06*natRNw^6 - 6.12758e-04*natRNw^5 + 2.41827e-02*natRNw^4 - 5.05612e-01*natRNw^3 + 5.91688e+00*natRNw^2 - 3.70536e+01*natRNw + 1.00265e+02);
                SwetW = 2*(Sref-(0.85*2*mac/(1+tratio))*fusd)*1.02;
                if SA >= 0 && SA < 10
                    Kw = 267.45*tcavg^4 - 70.235*tcavg^3 + 9.7108*tcavg^2 + 1.5765*tcavg + 1.0008;
                elseif SA >= 10 && SA < 15
                    Kw1 = 267.45*tcavg^4 - 70.235*tcavg^3 + 9.7108*tcavg^2 + 1.5765*tcavg + 1.0008;
                    Kw2 = 28.578*tcavg^3 - 1.6311*tcavg^2 + 1.8873*tcavg + 0.9989;
                    Kw = Kw1 - (Kw1-Kw2)/(10-15)*(10-SA);
                elseif SA >= 15 && SA < 20
                    Kw1 = 28.578*tcavg^3 - 1.6311*tcavg^2 + 1.8873*tcavg + 0.9989;
                    Kw2 = 25.476*tcavg^3 - 0.868*tcavg^2 + 1.7751*tcavg + 0.9995;
                    Kw = Kw1 - (Kw1-Kw2)/(15-20)*(15-SA);
                elseif SA >= 20 && SA < 25
                    Kw1 = 25.476*tcavg^3 - 0.868*tcavg^2 + 1.7751*tcavg + 0.9995;
                    Kw2 = 282.52*tcavg^4 - 70.398*tcavg^3 + 9.0241*tcavg^2 + 1.3817*tcavg + 0.9999;
                    Kw = Kw1 - (Kw1-Kw2)/(20-25)*(20-SA);
                elseif SA >= 25 && SA < 30
                    Kw1 = 282.52*tcavg^4 - 70.398*tcavg^3 + 9.0241*tcavg^2 + 1.3817*tcavg + 0.9999;
                    Kw2 = 162.49*tcavg^4 - 29.447*tcavg^3 + 4.508*tcavg^2 + 1.4476*tcavg + 1.0008;
                    Kw = Kw1 - (Kw1-Kw2)/(25-30)*(25-SA);
                elseif SA >= 30 && SA < 35
                    Kw1 = 162.49*tcavg^4 - 29.447*tcavg^3 + 4.508*tcavg^2 + 1.4476*tcavg + 1.0008;
                    Kw2 = 212.41*tcavg^4 - 48.243*tcavg^3 + 6.7731*tcavg^2 + 1.2485*tcavg + 1.0015;
                    Kw = Kw1 - (Kw1-Kw2)/(30-35)*(30-SA);
                elseif SA >= 35 && SA <= 40
                    Kw1 = 212.41*tcavg^4 - 48.243*tcavg^3 + 6.7731*tcavg^2 + 1.2485*tcavg + 1.0015;
                    Kw2 = 194.6*tcavg^4 - 38.877*tcavg^3 + 5.1437*tcavg^2 + 1.2233*tcavg + 1.0009;
                    Kw = Kw1 - (Kw1-Kw2)/(35-40)*(35-SA);
                end;
                fwing = Kw*Cfw*SwetW;
                ftail = (0.35+ 0.10/NoEng*Location)*fwing;
                
                % Fuselage
                RNf = RNperFT*fusl;
                natRNf = log(RNf);
                Cff = 1e-3*exp(6.42329e-06*natRNf^6 - 6.12758e-04*natRNf^5 + 2.41827e-02*natRNf^4 - 5.05612e-01*natRNf^3 + 5.91688e+00*natRNf^2 - 3.70536e+01*natRNf + 1.00265e+02);
                Swetf = 0.9*pi*fusl*fusd;
                LDratio = fusl/fusd;
                fig114 = csvread('Fig11.4data.txt');
                Kf = polyval(polyfit(fig114(:,1),fig114(:,2),3),LDratio);
                ffus = Kf*Cff*Swetf;
                
                % Nacelles/Pylon
                SwetN = 2.1*(Tengine)^.5*NoEng;
                fnac = 1.25*Cfw*SwetN;
                fpylon = 0.20*fnac;
                
                ftotal = (fwing+ftail+ffus+fnac+fpylon)*1.06;
                CDo = ftotal/Sref;
                e = 1/(1.035+0.38*CDo*pi*AR);
                fig6Ldg = csvread('fig6Landing.txt');
                CDoLanding = polyval(polyfit(fig6Ldg(:,2),fig6Ldg(:,1),3),CL/CLmaxldg);
                fig6To = csvread('fig6Takeoff.txt');
                CDoTakeoff = polyval(polyfit(fig6To(:,2),fig6To(:,1),3),CL/CLmaxto);
                
                %% Climb
                Vcl = 1.3*12.9/(ftotal*e)^.25*((WTO-.0175*WTO)/(0.5702*b))^.5;
                Mcl = Vcl/(0.592484*sqrt(1.4*1716.49*(-12.26+459.67)));
                TR = 0.5702*ftotal*Vcl^2/296+94.1/(0.5702*e*Vcl^2)*(WTO/b)^2;
                if strcmp(etype,'JT9D') == 1
                    TA1 = -23053*Mcl^5 + 73875*Mcl^4 - 91469*Mcl^3 + 60201*Mcl^2 - 30922*Mcl + 27143;
                    c1 = 7.604e-02*Mcl^2 + 3.471e-01*Mcl + 3.354e-01;
                    TA2 = 157192*Mcl^6 - 583328*Mcl^5 + 856431*Mcl^4 - 634808*Mcl^3 + 252828*Mcl^2 - 54627*Mcl + 20289;
                    c2 = 0.3698*Mcl + 0.3417;
                    cClimb = (c1+c2)/2;
                    TAcl = (TA1+TA2)/2;
                end;
                if strcmp(etype,'JT8D') == 1
                    TAcl = 36179*Mcl^6+122034*Mcl^5-150591*Mcl^4+81423*Mcl^3-14946*Mcl^2-2643.2*Mcl+7472;
                    cClimb = -0.9887*Mcl^6+2.1185*Mcl^5-0.7239*Mcl^4-1.216*Mcl^3+0.9944*Mcl^2+0.0909*Mcl+0.5905;
                end;
                if AdvancedEngineTec == 1
                    cClimb = cClimb*0.9;
                end;
                TA = Tengine/TSLT*TAcl;
                RClimb = 101*(NoEng*TA-TR)*Vcl/WTO;
                TimeCl = Altitude/RClimb; % Time to Climb (min)
                RangeCl = Vcl*TimeCl/60;
                WfCl = NoEng*TA*cClimb*TimeCl/60;
                
                %% Range
                W_0 = WTO-WfCl;
                W_1 = (1-Wf2Wto(z))*WTO;
                Cl_cruise = (W_0+W_1)/(2*Sref)/(1481*0.2360*MCruise^2);
                Cd = Cl_cruise^2/(pi*AR*e)+CDo+0.001;
                L2D = Cl_cruise/Cd;
                Trr = (W_0+W_1)/2/L2D;
                TrJT9D = Trr*TSLT/Tengine/NoEng;
                if strcmp(etype,'JT8D') == 1
                    c = 9.553e-15*TrJT9D^5-1.814e-10*TrJT9D^4+1.376e-6*TrJT9D^3-5.218e-3*TrJT9D^2+9.884*TrJT9D-7484;
                end;
                if strcmp(etype,'JT9D') == 1
                    c = 3.3144e-17*TrJT9D^4 - 1.3542e-12*TrJT9D^3 + 2.2789e-08*TrJT9D^2 - 1.6998e-04*TrJT9D + 1.0587; % Fit equation from 35000 ft chart
                end;
                if AdvancedEngineTec == 1
                    c = c*0.9;
                end;
                R_cruise = Vcruise/c*L2D*log(W_0/W_1);
                RangeTotal = RangeCl+R_cruise;
            end;
            if diff2 < tolRange || ClimbCheck == 1
                Wf2WtoConverge = 1;
                %% Check on Treq at Top of Climb
                CLIC = (W_0/Sref)/(1481*0.236*MCruise^2);
                CDinduced = CLIC^2/(pi*AR*e);
                CDdrag = CDo+CDinduced+0.0010;
                LLDD = CLIC/CDdrag;
                Treq = W_0/LLDD/NoEng;
                TrrJT9D = Treq*TSLT/Tengine;
                Tavil = 10296;
                if TrrJT9D > Tavil
                    WT = WT - tolWT;
                    ClimbCheck = 1;
                else
                    ClimbCheck = 0;
                    %% Climb Gradients
                    MaxCT = csvread('MaxCont T.txt');
                    fig6Ldg = csvread('fig6Landing.txt');
                    fig6To = csvread('fig6Takeoff.txt');
                    % 1st Segment
                    CLTO1 = CLmaxto/(1.2)^2;
                    CDo1 = polyval(polyfit(fig6To(:,2),fig6To(:,1),3),1/1.2^2);
                    CD1seg = CDo*2+CDo1+CLTO1^2/(pi*AR*e);
                    L2D1 = CLTO1/CD1seg;
                    Treq1 = WTO/L2D1;
                    V_1seg = 1.2*(296*WSto/0.953/CLTO1)^0.5;
                    M_1seg = V_1seg/(0.592484*sqrt(1.4*1718*(84+459.67)));
                    if strcmp(etype,'JT8D') == 1
                        Ta1 = Tengine/TSLT*(7606.2*M_1seg^2-9185.4*M_1seg+14401);
                    end;
                    if strcmp(etype,'JT9D') == 1
                        Ta1 = Tengine/TSLT*(37470*M_1seg^2 - 47380*M_1seg + 45368);
                    end;
                    Grad_1 = ((NoEng-1)*Ta1-Treq1)/WTO;
                    
                    % 2nd Segment
                    CD2seg = CDo+CDo1+CLTO1^2/(pi*AR*e);
                    L2D2 = CLTO1/CD2seg;
                    Treq2 = WTO/L2D2;
                    Grad_2 = ((NoEng-1)*Ta1-Treq2)/WTO;
                    
                    % 3rd Segment
                    if (0<=SA) && (SA<15)
                        figCleanWing1 = csvread('CLSA0.txt');
                        figCleanWing2 = csvread('CLSA15.txt');
                        CLmaxClean1 = interp1(figCleanWing1(:,1),figCleanWing1(:,2),tcavg,'spline');
                        CLmaxClean2 = interp1(figCleanWing2(:,1),figCleanWing2(:,2),tcavg,'spline');
                        CLmaxClean = CLmaxClean1+((SA/15)*(CLmaxClean1-CLmaxClean2));
                    end;
                    if (15<=SA) && (SA<35)
                        figCleanWing1 = csvread('CLSA15.txt');
                        figCleanWing2 = csvread('CLSA35.txt');
                        CLmaxClean1 = interp1(figCleanWing1(:,1),figCleanWing1(:,2),tcavg,'spline');
                        CLmaxClean2 = interp1(figCleanWing2(:,1),figCleanWing2(:,2),tcavg,'spline');
                        CLmaxClean = CLmaxClean1+((SA/20)*(CLmaxClean1-CLmaxClean2));
                    end;
                    if (SA>=35)
                        figCleanWing = csvread('CLSA15.txt');
                        CLmaxClean = interp1(figCleanWing(:,1),figCleanWing(:,2),tcavg,'spline');
                    end;
                    V_3seg = 1.2*(296*WSto/0.925/CLmaxClean)^.5;
                    M_3seg = V_3seg/(0.592484*sqrt(1.4*1716.49*(515.44)));
                    CL_3seg = CLmaxClean/(1.2)^2;
                    CD_3seg = CDo+CL_3seg^2/(pi*AR*e);
                    L2D3 = CL_3seg/CD_3seg;
                    Treq3 = WTO/L2D3;
                    if strcmp(etype,'JT8D') == 1
                        Ta3 = Tengine/TSLT*(-2172.3*M_3seg^4+120.42*M_3seg^3+8164.4*M_3seg^2-9004.9*M_3seg+12597);
                    end;
                    if strcmp(etype,'JT9D') == 1
                        Ta3 = Tengine/TSLT*(-2172.3*M_3seg^4+120.42*M_3seg^3+8164.4*M_3seg^2-9004.9*M_3seg+12597);
                    end;                    
                    Grad_3 = ((NoEng-1)*Ta3-Treq3)/WTO;
                    
                    % Approach
                    CLapp = CLmaxto/1.3^2;
                    CDoapp = polyval(polyfit(fig6To(:,2),fig6To(:,1),3),1/1.3^2);
                    CDapp = CDo+CDoapp+CLapp^2/(pi*AR*e);
                    L2Dapp = CLapp/CDapp;
                    WLDG = WSldg*Sref;
                    Treqapp = WLDG/L2Dapp;
                    Vappp = 1.3*(296*WSldg/0.953/CLapp)^.5;
                    Mappp = Vappp/(0.592484*sqrt(1.4*1716.49*(84+459.67)));
                    if strcmp(etype,'JT8D') == 1
                        Taapp = Tengine/TSLT*(-2172.3*Mappp^4+120.42*Mappp^3+8164.4*Mappp^2-9004.9*Mappp+12597);
                    end;
                    if strcmp(etype,'JT9D') == 1
                        Taapp = Tengine/TSLT*(-2172.3*Mappp^4+120.42*Mappp^3+8164.4*Mappp^2-9004.9*Mappp+12597);
                    end;
                    Grad_AP = ((NoEng-1)*Taapp - Treqapp)/WLDG;
                    
                    %Landing
                    CL_ldg = CLmaxldg/1.3^2;
                    CDoldg = polyval(polyfit(fig6Ldg(:,2),fig6Ldg(:,1),3),1/1.3^2);
                    CD_ldg = CDo+CDoldg+CDo+CL_ldg^2/(pi*AR*e);
                    L2Dldg = CL_ldg/CD_ldg;
                    Treqldg = WLDG/L2Dldg;
                    Mldg = Vapp/(0.592484*sqrt(1.4*1716.49*(84+459.67)));
                    if strcmp(etype,'JT8D') == 1
                        Taldg = Tengine/TSLT*(7606.2*Mldg^2-9185.4*Mldg+14401);
                    end;
                    if strcmp(etype,'JT9D') == 1
                        Taldg = Tengine/TSLT*(37470*Mldg^2 - 47380*Mldg + 45368);
                    end;
                    Grad_ldg = (NoEng*Taldg-Treqldg)/WLDG;
                    if NoEng == 2
                        if Grad_1 > 0.000 || Grad_2 > 0.024 || Grad_3 > 0.012 || Grad_AP > 0.021 || Grad_ldg > 0.032
                            start = 1;
                            GradConverge = 1;
                            break;
                        else
                            start = 0;
                            WT = WT-tolWT;
                            iteration = iteration + 1;
                        end;
                    elseif NoEng == 3
                        if Grad_1 > 0.003 || Grad_2 > 0.027 || Grad_3 > 0.015 || Grad_AP > 0.024 || Grad_ldg > 0.032
                            start = 1;
                            GradConverge = 1;
                            break;
                        else
                            start = 0;
                            WT = WT-tolWT;
                            iteration = iteration + 1;
                        end;
                    elseif NoEng == 4
                        if Grad_1 > 0.005 || Grad_2 > 0.030 || Grad_3 > 0.017 || Grad_AP > 0.0217|| Grad_ldg > 0.032
                            start = 1;
                            GradConverge = 1;
                            break;
                        else
                            start = 0;
                            WT = WT-tolWT;
                            iteration = iteration + 1;
                        end;
                    else
                        start = 1;
                    end;
                end;
            else
                start = 1;
                Wf2WtoConverge = 0;
            end;
        end;
        if GradConverge == 1
            break;
        end;
    end;
    if fail ~= 1
        %% Direct Operating Cost
        D1 = Range*1.15;
        TGM = D1/(11.866+0.040669*D1)/60;
        TCL = TimeCl/60;
        TAM = .10;
        TCR = (D1*1.02+20-RangeCl*1.15)/(1.15*Vcruise);
        TB = TGM+TCL+TCR+TAM;
        VB = D1/TB;
        FB = WfCl+Trr*c*(TCR+TAM);
        
        %% Flying Operating Cost
        % Flight Crew
        P = WPL/2000;
        perBLKhr = 17.849*(1.15*Vcruise*WTOal/10^5)^0.3+40.83;
        CTM_fc = perBLKhr/(VB*P);
        % Fuel/Oil
        CFT = FuelCost/6.4;
        COT = 2.15;
        CTM_fo = (1.02*FB*CFT+NoEng*COT*TB*.135)/(D1*P);
        % Hull Insurance
        Wa = WTOal*(1-Wf2Wto(length(Wf2Wto)))-WPL-WPP*WTOal;
        CA = 2.4e6+87.5*Wa;
        if AdvancedEngineTec == 1
            supppp = 1.1;
        else
            supppp = 1;
        end;
        CE = (590000+Tengine*16)*supppp;
        CT = CA+NoEng*CE;
        IRa = .01;
        U = 630+4000/(1+1/(TB+0.5));
        CTM_hi = IRa*CT/(U*VB*P);
        % DM Labor-Airplane
        KFHA = 4.9169*log10(Wa/10^3)-6.425;
        KFCA = 0.21256*(log10(Wa/10^3))^3.7375;
        TF = TB-TGM;
        RL = 8.60;
        CTM_al = (KFHA*TF+KFCA)/(VB*TB*P)*RL*(1+0.29*(MCruise-1))^1.5;
        % DM Airframe Material
        CFHA = 1.5994*CA/10^6+3.4263;
        CFCA = 1.9229*CA/10^6+2.2504;
        CTM_am = (CFHA*TF+CFCA)/(VB*TB*P);
        % DM Engine-Labor
        KFHE = NoEng*(Tengine/10^3)/(0.82715*(Tengine/10^3)+13.639);
        KFCE = .20*NoEng;
        CTM_el = (KFHE*TF+KFCE)/(VB*TB*P)*RL;
        % DM Engine-Material
        CFHE = (28.2353*CE/10^6-6.5176)*NoEng;
        CFCE = (3.6698*CE/10^6+1.3685)*NoEng;
        CTM_em = (CFHE*TF+CFCE)/(VB*TB*P);
        % DM Total Maintenance-Burdened
        CTM_to = (CTM_al+CTM_am+CTM_el+CTM_em)*2;
        % Depreciation
        CTM_de = 1/(VB*P)*((CT+0.06*(CT-NoEng*CE)+0.3*NoEng*CE)/(14*U));
        TONMILE = CTM_fc+CTM_fo+CTM_hi+CTM_to+CTM_de;
        TotalDOC = [CTM_fc CTM_fc/TONMILE; CTM_fo CTM_fo/TONMILE;
            CTM_hi CTM_hi/TONMILE; CTM_to CTM_to/TONMILE; CTM_de CTM_de/TONMILE];
        PAXMILE = TONMILE*P/NoPAX;
        sizeWf22 = size(Wf2Wto);
        Wf2Wto = Wf2Wto(sizeWf22(1,2));
        fprintf('\nFor a SA = %g AR = %g and Wing Type = %s',SA,AR,wtype);
        fprintf('\nConverged CL = %g', CL);
        fprintf('\nConverged Wf2Wto = %g', Wf2Wto);
        fprintf('\nConverged WT = %g', WT);
        fprintf('\nTotal Flat Plate Area = %g', ftotal);
        fprintf('\nSref = %g   b = %g   m.a.c. = %g',Sref,b,mac);
        fprintf('\nRN_wing = %g  Cf_wing = %g',RNw,Cfw);
        fprintf('\nRN_fuse = %g  Cf_fuse = %g',RNf,Cff);
        fprintf('\nTotal DOC per PAXMILE = %g', PAXMILE);
        fprintf('\nTotal DOC per TONMILE = %g\n', TONMILE);
        a(q,:) = [SA AR];
        aa(q,:) = [TONMILE PAXMILE];
        Store(q,:) = [NoPAX SA AR tratio {wtype} NoEng Location NoAbreast NoAisles Wf2Wto CL CLmaxto CLmaxldg WT WTO Sref b tcavg ftotal e L2D TONMILE PAXMILE];
        Fuell(q,:) = [FuelCost TONMILE];
    end;
    if strcmp(VaryAR,'yes')
        AR = AR+Increment;
    elseif strcmp(VarySA,'yes')
        SA = SA+Increment;
    elseif strcmp(VaryFC,'yes')
        FuelCost = FuelCost+Increment;
    end;
end;

%% Output to Excel
if output == 1
    if strcmp(VaryFC,'yes')
        xlswrite(filename,Fuell,sheet,strcat('A',num2str(Row)));
    else
        if addTitle == 1
            Title = [{'NoPAX'} {'SA'} {'AR'} {'tratio'}	{'wtype'} {'NoEng'}	{'Location'} {'NoAB'} {'NoAisle'} {'Wf/Wto'} {'CL'}	{'CL_TO'} {'CL_LDG'} {'W/T'} {'W_TO'} {'Sref'} {'b'} {'t/c'} {'f'} {'e'} {'L/D'} {'TONmile'} {'PAXmile'}];
            xlswrite(filename,Title,sheet,strcat('A','1'));
        end;
        xlswrite(filename,Store,sheet,strcat('A',num2str(Row)));
        [rows,col] = size(Store);
        endRow = rows+Row;
        fprintf('Outputted values to Row %g to %g.',Row,endRow);
        if addTitle ==1
            fprintf(' Added a title as well');
        end;
        fprintf('\n');
    end;
end;
%plot(a(:,1),aa(:,1))
toc