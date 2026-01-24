%% Problem 1 – Dynamic Programming
clear; clc;

% N = 20 stages  → k = 0,...,19
% Feasible states x in {1,...,c}
% Feasible controls u in {1,...,d}
% Random variable w in {1,...,6} with equal probability

% Parameters:
N = 20;
c = 12;
d = 8;
x0 = 2;

states   = 1:c;
controls = 1:d;
wvals    = 1:6;
pw       = 1/6;

% Terminal cost:
% g_N(x_N) = 3*x_N + min(x_N, 3)
J = zeros(N+1, c);
for x = states
    J(N+1, x) = 3*x + min(x, 3);
end

% BACKWARD DP RECURSION based on Bellmans equation:
% J*_k(x) = min_u  E_w[ g(x,u,w) + J*_{k+1}( f(x,u,w) ) ]
% Where:
%   System equation f(x,u,w):
%   x_{k+1} = min{ max[x + 2u − 2w − 3, 1], c }
%   Stage cost:
%   g(x,u,w) = x + 2u − 5 * min(x,5)
for k = N:-1:1      % iterate starts from 0
    for x = states
        
        bestVal = inf;   % declaracion
        
        for u = controls
            % E_w[ g(x,u,w) + J_{k+1}(x_next) ]
            Ek = 0;
            
            for w = wvals
                % System equation:
                % x_{k+1} = min{ max[x + 2u − 2w − 3, 1], c }
                xn = x + 2*u - 2*w - 3;
                xn = max(1, min(c, xn));
                
                % Stage cost:
                % g(x,u,w) = x + 2u − 5 * min(x,5)
                g = x + 2*u - 5*min(x, 5);
                
                % Expected value accumulation
                Ek = Ek + pw * (g + J(k+1, xn));
            end
            
            % Minimization over controls u
            bestVal = min(bestVal, Ek);
        end
        
        % Store optimal value J_k(x)
        J(k, x) = bestVal;
    end
end

% Result for J(x0) - x0 = 2, c = 12, d = 8
fprintf("J(x0 = %d) = %.6f\n", x0, J(1, x0));
