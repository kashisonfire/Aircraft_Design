%% Tail Sizing
% Run Main.m, Run CoG.m, and input a Wing before running this code
clear Table im E
Wing =66.2;
CoG2 = 82.5;
E = 0:1:5;
Names = {'E' 'LengthVT' 'L_VT' 'CR_VT' 'CT_VT' 'mac_VT' 'CR_HT' 'CT_HT' 'mac_HT'};
for i = 1:length(E)
    %% Vertical Tail
    L_VT = fusl+E(i)-Wing-y_mac;
    S_VT = 0.08*Sref*b/L_VT;
    AR_VT = 2.5;
    b_VT = sqrt(AR_VT*S_VT);
    CR_VT = 2*S_VT/(b_VT*(1+tratio));
    CT_VT = tratio*CR_VT;
    mac_VT = 2/3*(CR_VT+CT_VT-CR_VT*CT_VT/(CR_VT+CT_VT));
    
    %% Horizontal Tail
    L_HT = L_VT;
    S_HT = 1.1*Sref*mac/L_HT;
    AR_HT = 5;
    b_HT = sqrt(AR_HT*S_HT);
    CR_HT = 2*S_HT/(b_HT*(1+tratio));
    CT_HT = tratio*CR_HT;
    mac_HT = 2/3*(CR_HT+CT_HT-CR_HT*CT_HT/(CR_HT+CT_HT));
    pos = E(i);
    LengthVT(i) = CoG2+L_VT;
    im(i,:) = [pos LengthVT(i) L_VT CR_VT CT_VT mac_VT CR_HT CT_HT mac_HT];
end;
Table = [Names;num2cell(im)];