function [dydt] = p09_sys(t,y,a,b)
% ODE system 
% y1 = x;
% y2 = p;
dydt = [a*y(1)-b^2*y(2);
        -a*y(2)];
end

