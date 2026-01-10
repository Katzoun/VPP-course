function [res] = p09_cond(ya,yb,x0,q)
% ya - vector of left conditions (t=0)
% yb - vector of right conditions (t=T)
res = [ya(1)-x0;
       yb(2)-q*yb(1)];
end

