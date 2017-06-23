function [z,exitflag,iter,lam,t] = quad_wright(H,f,A,b,maxiter,tol,verbose,z0,lam0,t0)
% Solve quadratic programming problem using Wright's (1997) Method
% Minimise 1/2x'Hx + f'x
% Subject to: Ax <= b 
%
% [z,exitflag,iter,lam,t] = quad_wright(H,f,A,b,maxiter,tol,verbose,z0,lam0,t0)

% Reference: S. J. Wright, "Applying New Optimization Algorithms to Model
% Predictive Control," in Chemical Process Control-V, CACHE, AIChE 
% Symposium, 1997, pp. 147-155.


% Length of constraint matrix
mc = length(b); 
% Number of decision variables
ndec = length(f);

% Default Args
if(nargin < 10), t = []; else t = t0; end
if(nargin < 9), lam = []; else lam = lam0; end
if(nargin < 8), z = []; else z = z0; end; 
if(nargin < 7 || isempty(verbose)), verbose = 0; end
if(nargin < 6 || isempty(tol)), tol = 1e-6; end
if(nargin < 5 || isempty(maxiter)), maxiter = 200; end
% Test for Warm Start
pmax = max(max(abs([H f; A b])));
if(pmax > 1)
    WARMVAL = sqrt(pmax);
else
    WARMVAL = 0.5;
end
if(isempty(z)) % cold
    z = zeros(ndec,1); 
    lam = WARMVAL*ones(mc,1); 
    t = WARMVAL*ones(mc,1);   
    wmode = 0;
elseif(isempty(lam)) % just primal
    lam = WARMVAL*ones(mc,1); 
    t = WARMVAL*ones(mc,1);
    wmode = 1;
elseif(isempty(t)) % primal + dual
    t = WARMVAL*ones(mc,1); 
    wmode = 2;
else % all
    wmode = 3;
end

% Default Values
sigma = 0.1; ascale = 1; inftol = tol*10;
exitflag = 1; cholfail = 0;
mu = t'*lam / mc; At = A';
mr2_1 = 100; mr2 = 10;
% Initial Residuals
r1 = -H*z - At*lam - f;
r2 = -A*z + b;
% Linsolve options
opU.UT = true; opUT.UT = true; opUT.TRANSA = true;

if(verbose)
    fprintf('-----------------------------------------------\n');
    fprintf('QuadWright QP Solver [MATLAB Double Version]\n');
    switch(wmode)
        case 3, fprintf(' Warm Start: Primal + Dual + Slack\n');
        case 2, fprintf(' Warm Start: Primal + Dual\n');
        case 1, fprintf(' Warm Start: Primal\n');
    end
    fprintf(' %4d Decision Variables\n %4d Inequality Constraints\n',ndec,mc);
    fprintf('-----------------------------------------------\n');
    fprintf('iter           phi             mu        sigma        alpha       max(r1)       max(r2)\n');
end 

% Begin Searching
for iter = 1:maxiter
    % Create common matrices    
    ilam = 1./lam;
    ilamt = ilam.*t;
    lamt = lam./t;
    mesil = mu*sigma.*ilam;    
    IGA = bsxfun(@times,A,lamt); % matrix * diagonal matrix
    igr2 = lamt.*(r2 - mesil);

    % Solve
    [R,p] = chol(H+At*IGA);
    if(~p)
        del_z = linsolve (R, linsolve (R, (r1+At*igr2), opUT), opU); % exploit matrix properties for solving
    else % Not Positive Definite 
        if(verbose), fprintf(2,'\b (Cholesky Failed)\n'); end
        del_z = (H+At*IGA)\(r1+At*igr2);
        cholfail = cholfail + 1;        
        if(cholfail > 2)
            exitflag = -2;
            break; 
        end
		% Pull back maximum step size
        ascale = ascale - 0.1;
    end
    del_lam = -igr2 + IGA*del_z;
    del_t = -t + mesil - ilamt.*del_lam;
    
    % Decide on suitable affine step-length
    duals = [lam;t];
    delta = [del_lam;del_t];
    index = delta < 0; 
    if(any(index))
        alpha = 0.9995*min(-duals(index)./delta(index)); % solves for min ratio (max value of alpha allowed)
    else
        alpha = 0.999995;
    end
    % Check for numerical problems (alpha will go negative)
    if(alpha < eps(1)), exitflag = -3; break; end
    % Local Scaling
    alpha = alpha*ascale;
    % Increment
    lam = lam + alpha*del_lam;
    t = t + alpha*del_t;
    z = z + alpha*del_z;

    % Update residuals
    r1 = (1-alpha)*r1; %equiv to r1 = -H*z - At*lam - f
    r2 = -A*z + b;
    % Complementarity Gap
    muold = mu;
    mu = t'*lam / mc;
    % Infeasibility Phi
    mr2_2 = mr2_1; mr2_1 = mr2;
    mr1 = max(abs(r1)); mr2 = max(abs(r2-t));
    phi = (max([mr1,mr2]) + t'*lam)/pmax;
    if(verbose)
        fprintf('%3d  %13.5g  %13.5g  %11.5g  %11.5g   %11.5g   %11.5g\n',iter,phi,mu,sigma,alpha,mr1,mr2);
    end
    % Check for NaNs
    if(isnan(mu) || isnan(phi)), exitflag = -3; break; end
    % Check Convergence
    if(mu <= tol && phi<=tol)
        exitflag = 1;
        if(verbose)
            fprintf('-----------------------------------------------\n');
            fprintf(' Successfully solved in %d Iterations\n Final phi: %d, mu %g [tol %g]\n',iter,phi,mu,tol);
            fprintf('-----------------------------------------------\n');
        end
        return
    end
    % Check For Primal Infeasible
    if(iter > 6 && mr2/pmax > tol/10)
        if(abs(mr2-mr2_1)/mr2 < inftol && abs(mr2_1-mr2_2)/mr2_1 < inftol)
            if(verbose), fprintf(2,'\b (Primal Infeasibility Detected)\n'); end  
            exitflag = -4;
            break;
        end  
    end  
    % Solve centering parameter (Mehrotra's Heuristic)
    sigma = min((mu/muold)^3,0.99999);
end

% If here, either bailed on Cholesky or Iterations Expired
if(exitflag==1), exitflag = -1; end %expired iters    
% Optional Output
if(verbose)
    fprintf('-----------------------------------------------\n');
    switch(exitflag) 
        case -1, fprintf(' Maximum Iterations Exceeded\n');
        case -2, fprintf(2,' Failed - Cholesky Factorization Reported Errors\n');
        case -3, fprintf(2,' Failed - Numerical Errors Detected\n'); 
        case -4, fprintf(2,' Failed - Problem Looks Infeasible\n');
    end
    fprintf(' Final phi: %g, mu %g [tol %g]\n',phi,mu,tol);
    fprintf('-----------------------------------------------\n');
end

