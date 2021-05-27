function  [ModeShapes,Slopes] = ReadBModOut(InputFile,MaxMode,MSIDs)
% function: read BModes output results
% Input
% InputFile: Output file of BModes (.out file), string
% MaxMode:   Number of modes in the .out file, integer
% MSIDs  : IDs of the mode shapes (1st & 2nd FA, 1st & 2nd SS), Vector

% Output
% ModeShapes: Shapes of the 4 modes, xx*5 matrix
% Slopes: 1*4 vector, slopes of each mode
%
% Written by: Yang Yang in Ningbo University
% Date:     : 25th May 2021
% email: yangyang1@nbu.edu.cn

frd = fopen(InputFile,'rt');
CurrentMode = 0;


while 1
    str = fgetl(frd);
    TempLen = length(str);
    if TempLen>0 & str(1:8) == 'span_loc'
       CurrentMode = CurrentMode +1;
       fgetl(frd);
       j=1;
       while 2
             str =  fgetl(frd);
             TempData = sscanf(str,'%f',5); 
             Ht(j) = TempData(1);
             MSDisp_ss(j,CurrentMode) = TempData(2);
             MSSlop_ss(j,CurrentMode) = TempData(3);
             MSDisp_fa(j,CurrentMode) = TempData(4);
             MSSlop_fa(j,CurrentMode) = TempData(5);
            if Ht(j) == 1.0
                break 
            end
            j=j+1;
       end
       if CurrentMode == MaxMode
          break 
       end
    end
    
end
ModeShapes(:,1) =  Ht';
for i=1:2
    ModeShapes(:,i+1) = MSDisp_fa(:,MSIDs(i));
    Slopes(i) = MSSlop_fa(1,MSIDs(i));
end
for i=1:2
    ModeShapes(:,i+3) = MSDisp_ss(:,MSIDs(i+2));
    Slopes(i+2) = MSSlop_ss(1,MSIDs(i+2));
end


        

