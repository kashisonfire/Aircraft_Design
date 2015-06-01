clc;
close all;
filename = 'Output.xlsx';
sheet = '1';
AspectRatio = xlsread(filename,sheet,'C2:C23');
DOC_TONmile = xlsread(filename,sheet,'V2:V23');
data = [AspectRatio DOC_TONmile];
data = sortrows(data,1);
ARtest = 7.75:0.001:8.35;
DOC = interp1(data(:,1),data(:,2),ARtest,'cubic');
plot(ARtest,DOC,AspectRatio,DOC_TONmile,'x')
for i=2:length(ARtest)
    if DOC(i) < DOC(i-1)
        OptDOC = DOC(i);
        OptAR = ARtest(i);
    end;
end;
disp(OptAR)
disp(OptDOC)