clear; clc
alpha = 0.8; x0 = 1;

u0 = x0 / (1 + alpha^2 + alpha^4);
x1 = x0 - u0;
u1 = x1 / (1 + alpha^2);
x2 = x1 - u1;       % = u2
u2 = x2;

Jstar = -2*( sqrt(u0) + alpha*sqrt(u1) + alpha^2*sqrt(u2) );

[u0, u1, u2, Jstar]