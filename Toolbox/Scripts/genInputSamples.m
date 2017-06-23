function [T_,U_,T,U] = genInputSamples(N,Ts,Uss,varargin)

% Correct Uss dimensions to nux1
if size(Uss,1) == 1
    Uss = Uss';
end
Uss = Uss(:,1);

% Define ranges
if nargin > 3 % Maximum samples specified by user
    n_max = varargin{1};
else
    n_max = 20; % Default of 20 maximum sample time (generally it is tp/Ts, where tp is process time, Ts is sampling time)
end

if nargin > 4 % Max. input change is specified by user
    dU_max = varargin{2}/2;
else
    dU_max = 0.1*Uss/2; % Default of +/- 5% steady-state value
end

if nargin > 5 % fbacktrace specified
    fbacktrace = varargin{3};
else
    fbacktrace = true; % Default backtrace as true (2N samples will be generated)
end

nu = length(Uss);

%% Generate LHS points
% Note the LHS generated will be used as follows
% First column translates into dt
% Second column translates into u_idx (which input to vary)
% Third column translates into du (corresponding to u_idx)

X = getLHSSamples(N,3);

T = zeros(1,N+1); U = zeros(nu,N+1);
U(:,1) = Uss;
n_bias = 1*nu; % Denotes how many steps to keep at SS initially
if N <= n_bias
    n_bias = nu;
end
dt_bias = Ts; % Denotes time each-step will be held at initially

for k = 1:N
    %dt = ceil(n_max*X(k,1))*Ts;
    dt = n_max*Ts/ceil(3*X(k,1)); % Divide into 3 random intervals
    %dt = n_max*Ts; % Uncomment this if all excitation interval is desired to be same 
    if dt <= 0 % very rare occasion when rand() yields 0.0000....
        dt = Ts;
        disp('dt = 0 encountered. This is a rare event. Forcing this to be Ts...');
    end
    if k <= n_bias
        dt = dt_bias;
    end
    T(k+1) = T(k) + dt;
    
    u_idx = ceil(nu*X(k,2));
    if u_idx <= 0 % very rare occasion when rand() yields 0.0000....
        u_idx = 1;
        disp('u_idx = 0 encountered. This is a rare event. Forcing this to be [1]...');
    end
    %u_idx = rem(k-1,nu) + 1; % Rotate indices e.g. 1,2,3,1,2,3,1,2,3...etc.
    
    if k > 1
        U(:,k) = U(:,k-1);
    end
    if k > n_bias
        U(u_idx,k) = Uss(u_idx) + dU_max(u_idx)*(2*X(k,3)-1);
    end
end
U(:,k+1) = U(:,k);

if fbacktrace
    U(:,end) = [];
    U = [U U(:,end-1:-1:1)]; T = [T zeros(1,N)];
    for k = 1:N
        T(N+1+k) = T(N+k) + (T(N+2-k) - T(N+1-k));
    end
    U = [U U(:,end) U(:,end)];
end

%% Change everything to fixed discrete intervals
T_ = 0:Ts:T(end);
U_ = zeros(nu,length(T_));
n = 1;
for k = 1:length(T_)
    if T_(k) >= T(n)
    %if abs(T_(k)-T(n)) < 1e-5 % if T_(k) == T(n)
        U_(:,k) = U(:,n);
        n = n + 1;
    else
        U_(:,k) = U_(:,k-1);
    end
end
