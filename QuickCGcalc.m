wing_aft = 48.9;
wing_forward = 47.25;
wing_mac25 = 48.075;

CG_forward = interp1(COG_npf(:,1),mac25,wing_forward,'linear');
CG_aft = interp1(COG_npf(:,1),mac25,wing_aft,'linear');
CG_mac25 = interp1(COG_npf(:,1),mac25,wing_mac25,'linear');

nose2main = 54; %approx. of distance
