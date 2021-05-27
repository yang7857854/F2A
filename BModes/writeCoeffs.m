function [] = writeCoeffs(Coef,fileNameOut)
% fucntion: write mode shape coefficients
%Input
% Coef: coefficients obtained from "getModeShapeCoef"
% fileNameOut: output file name
%
% Written by: Yang Yang in Ningbo University
% Date:     : 25th May 2021
% email: yangyang1@nbu.edu.cn

fid1 = fopen(fileNameOut,'wt');
    for i=1:5
        TempS = num2str(i+1);
        fprintf(fid1,'%13.7f \t %s \n',Coef(i), strcat('TwFAM1Sh(',TempS,') - Mode 1, coefficient of x^',TempS,' term'));
    end
    for i=1:5
        TempS = num2str(i+1);
        fprintf(fid1,'%13.7f \t %s \n',Coef(i+5), strcat('TwFAM2Sh(',TempS,') - Mode 2, coefficient of x^',TempS,' term'));
    end
    fprintf(fid1,'%s \n','---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------');
    for i=1:5
        TempS = num2str(i+1);
        fprintf(fid1,'%13.7f \t %s \n',Coef(i+10), strcat('TwSSM1Sh(',TempS,') - Mode 1, coefficient of x^',TempS,' term'));
    end
    for i=1:5
        TempS = num2str(i+1);
        fprintf(fid1,'%13.7f \t %s \n',Coef(i+15), strcat('TwSSM2Sh(',TempS,') - Mode 2, coefficient of x^',TempS,' term'));
    end
fclose('all');  
