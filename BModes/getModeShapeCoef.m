 function Coef = getModeShapeCoef(ModeShapes,Slopes,MethodID)
% function: calcuate mode shape polyfitting coefficients
% Input
% ModeShapes: Shapes of the 4 modes, xx*5 matrix
% Slopes: 1*4 vector, slopes of each mode
% MethodID : 1 for "Projection", 2 for "ImprovedDirect" 
%
% Output
% Coeffs: Coefficients of the given modes

% Written by: Yang Yang in Ningbo University
% Date:     : 25th May 2021
% email: yangyang1@nbu.edu.cn
%MethodID = 2; % 1 for "Projection", 2 for "ImprovedDirect" 
x = ModeShapes(:,1);
NumModes = length(ModeShapes(1,:))-1;
Coef = zeros(5*NumModes,1);
for i=1:NumModes
    x = x-x(1);
    y=ModeShapes(:,i+1);
    % projection method
    f1 = 1.0e-4*abs(x(2)-x(1))/(y(2)-y(1));
    f2 = Slopes(i);
    factor_y = y.*f1;
    project_x1 = x.*cos(atan(f1*f2))+(factor_y-factor_y(1)).*(sin(atan(f1*f2)));
    Norm_x1 = project_x1./(max(project_x1));
    project_y1 = -x.*sin(atan(f1*f2))+(factor_y-factor_y(1)).*(cos(atan(f1*f2)));
    Norm_y1 = project_y1./(project_y1(end));
    x_proj_s2 = Norm_x1.^2;
    x_proj_s3 = Norm_x1.^3;
    x_proj_s4 = Norm_x1.^4;
    x_proj_s5 = Norm_x1.^5;
    x_proj_s6 = Norm_x1.^6;
    % improved directed method
    y_direct = y-y(1)-x.*f2;
    Norm_x = x./(max(x));
    x_ipdr_s2 =Norm_x.^2; 
    x_ipdr_s3 =Norm_x.^3;
    x_ipdr_s4 =Norm_x.^4;
    x_ipdr_s5 =Norm_x.^5;
    x_ipdr_s6 =Norm_x.^6;
    if  (MethodID == 1)
        X = [x_proj_s2 x_proj_s3 x_proj_s4 x_proj_s5 x_proj_s6];
        Y = project_y1;
    elseif(MethodID == 2)
        X = [x_ipdr_s2 x_ipdr_s3 x_ipdr_s4 x_ipdr_s5 x_ipdr_s6];
        Y = y_direct;
    end
    
    b = regress(Y,X);
    b = b./(sum(b));
    b(5) = 1.00000000-sum(b(1:4));
    Coef((i-1)*5+1:(i-1)*5+5)=b;
    project_y_polyfit = b(1)*x_proj_s2+b(2)*x_proj_s3+b(3)*x_proj_s4+b(4)*x_proj_s5+b(5)*x_proj_s6;
    figure;
    xlim([0,1]);ylim([0,1])
    plot(Norm_x1,Norm_y1,'bv');
    hold on;
    if  (MethodID == 1)
        plot(Norm_x1,project_y_polyfit,'r*')
        legend('Orignal data','Projection data','Location','NorthWest')
    elseif(MethodID == 2)
        plot(Norm_x1,project_y_polyfit,'r*')
        legend('Orignal data','Improved Direct data','Location','NorthWest')
    end
end


