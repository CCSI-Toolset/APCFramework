function [x,exitflag,iter] = quad_hildreth(H,f,A,b,maxiter,tol,z0)
% Solve quadratic programming problem using Hildreth's Method
% Minimise J(x) = 1/2x'Hx + f'x
% Subject to: Ax <= b 

% Based on parts by Liuping Wang from the book "Model Predictive Control
% System Design and Implementation Using MATLAB"

% Reference: D. A. Wismer & R. Chattergy "Introduction to Nonlinear 
% Optimization, a Problem Solving Approach," North-Holland, New York,
% 1978

iter = 0; exitflag = 1; % assume will end normally  

% Cholesky Factorization of H
R = chol(H);

% Global Minimum
x = -R\(R'\f); % -H\f;  

% Otherwise solve using dual mthod
P = A*(R\(R'\A')); % A/H * A';
K = b - A*x;
m = length(K); % # of constraints: should be a vector

% Initial condition
lambda = zeros(size(b)); % 
lambdap = lambda+tol*100; 

% Solve element by element
while (norm(lambda-lambdap) > tol) 
    iter = iter + 1; 
    lambdap = lambda;
    for i = 1:m
        la = -(P(i,:)*lambda - P(i,i)*lambda(i,1) + K(i,1))/P(i,i);
        lambda(i,1) = max(0,la);
    end

    if(iter > maxiter) % Solution timed out  
        exitflag = -1;
        break;
    end
end

x = x - (R\(R'\A'))*lambda; %  % x - H\A'*lambda;
end
