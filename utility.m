function [ utils ] = utility (cons,L)

%This function takes consumption as an argument and returns utility. The
%utility functin is CRRA except that we add a very small number (eps) to
%consumption so that the computer can deal wi

global gamma nu hrsWrk

if cons<=0 
   error('Error in utility. Consumption is <=0');
end
                                
l=(L)*(hrsWrk)+(~L);
if gamma == 1
    utils = log(cons^nu*l^(1-nu));
else
    utils = ((cons^nu*l^(1-nu))^(1-gamma)  )/(1-gamma);
end

end

