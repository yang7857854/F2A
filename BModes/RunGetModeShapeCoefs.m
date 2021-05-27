clear all
close all
BModesOut = 'OC3Hywind.out';
MSIDs = [8 11 7 10]; % corresponding to 1st FA, 2nd FA, 1st SS and 2nd SS, respectively.
MaxMode = 20;        % Number of mode in the .out file
MethodID = 1;        % 1 for "Projection" method, 2 for "Improved Direct" method
% get the mode shapes and slope at the bottom
[ModeShapes,Slopes] = ReadBModOut(BModesOut,MaxMode,MSIDs);
% calculate the shape coefficients
Coeffs = getModeShapeCoef(ModeShapes,Slopes,MethodID);
% output the results
fileNameOut = strcat(BModesOut(1:length(BModesOut)-4),'_ModeShapeCoefs.txt');
writeCoeffs(Coeffs,fileNameOut)
disp(strcat('Check mode shape coefficients in: ',fileNameOut))