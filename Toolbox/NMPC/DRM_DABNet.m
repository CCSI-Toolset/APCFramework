classdef DRM_DABNet < handle
% Class definition for DABNet based dynamic reduced models
%
% Written and developed by:
% Priyadarshi Mahapatra, URS Corporation / National Energy Technology Laboratory
% Jinliang Ma, URS Corporation / National Energy Technology Laboratory
%
% For details of D-RM Builer application and D-RMs, refer to CCSI's D-RM Builer user's manual
% 
% For detailed description of the mathematical formulation and concepts, please refer to:
% G.B. Sentoni and L.T. Biegler, State-Space Nonlinear Process Modeling: Identification and Universality, AIChE Journal, Vol. 44, No. 10, October 1998
%
% Revision List
% 2013/10 - Initial version containing constructor, initialize and
%           evalNextStep functions. NN method defined within the class.
% 2013/11 - [Major Revision] Inserted prediction capability for APC through
%           predictNextStep function
% 2014/03 - NN function now part of DRM_ANN class
% 2014/07 - Added additional functions for UKF calculation (external to this class)
% 2014/08 - Modified evalNextStep and predictNextStep functions to accept M control moves
% 2014/10 - [Major Revision] Added derivative formulation
% 2015/03 - Code revision for computational speed including
%           a. Avoid fixed derivatives (dxannP_dduM) evaluations everytime
%              predictNextStep function is invoked except when P, M changes
%              by introducing derivativeP & derivativeM properties
%           b. Preallocation of memory especially 3D-matrix yiP_duM
% 2015/05 - Added A_ann, B_ann as properties which are evaluated within constructor

    properties
        dt;             % sampling time interval
        Ts;             % alias of sampling time interval
        u_name;         % names of all input variables including those not varied, cell array of strings
        y_name;         % names of all output variables including those not varied, cell array of strings
        u_unit;         % units of input variables
        y_unit;         % units of output variables
        t_unit;         % unit of time
        u_idx;          % indices of varied inputs among "plant" inputs, integer array
        y_idx;          % indices of varied outputs among "plant" outputs, integer array
        nu;             % number of input variables
        ny;             % number of output variables
        nx_ann;         % vector of number of states for each output [ny]
        u_mean;         % vector of means of input variables, size [nu]
        y_mean;         % vector of means of output varaibles, size [ny]
        u_sigma;        % vector of standard deviations of input variables, size [nu]
        y_sigma;        % vector of standard deviations of output varaibles, size [ny]
        A;              % 2-D cell array of A matrices, size {ny,nu}
        B;              % 2-D cell array of B vectors, size {ny,nu}
        ANN;            % 1-D array of artificial neural networks, size [ny]
        x;              % 2-D cell array of decoupled state-space vector, size {ny,nu}
        y;              % output vector, size [ny]
        dy_du;          % jacobian matrix, size [ny,nu]
		A_ann;          % 1-D cell array of augmented A matrices, size {ny}
		B_ann;          % 1-D cell array of augmented B matrices, size {ny}
        dxannP_dduM;    % array of ny cells, containing dxann/du information for each output
        derivativeP = 0; derivativeM = 0;
        % Variables related to state-estimation
        fStateEstimation = false;
        fIncludeStates = false;
        d_idx;          % indices of estimated disturbances among model inputs, integer array;
        nd;
        d;              % estimated disturbances
        Mk; Pk; SDk;
        KF;
        FilterMethod = 'EKF'; % Choose between 'UKF' or 'EKF'
        SEFormulation = 'Random-Step Input Disturbance';
        qd; qx; r;
        Qd; Qx;
        Q; R;
        % For debugging
        %jacobian_num; jacobian_ana;
        %hessian_num; hessian_ana;
    end

    methods
        function obj = DRM_DABNet(Ts,A,B,NN,u_mean,y_mean,u_sigma,y_sigma,u_idx,y_idx,u_name,y_name,u_unit,y_unit,t_unit)
            % constructor, parameters of the function are assigned by a MATLAB script exported by D-RM Builder
            obj.dt = Ts; obj.Ts = Ts;
            obj.A = A; obj.B = B;
            obj.nu = size(obj.A,2); obj.ny = size(obj.A,1);
            obj.ANN = DRM_ANN.empty(obj.ny,0);
            for i=1:obj.ny
                obj.ANN(i) = DRM_ANN(NN(i));
            end
            obj.u_mean = u_mean'; obj.y_mean = y_mean';
            obj.u_sigma = u_sigma'; obj.y_sigma = y_sigma';
            % default I/O index: u_idx = [1 2 3 ... nu]; y_idx = [1 2 3 ... ny];
            obj.u_idx = 1:obj.nu; obj.y_idx = 1:obj.ny;
            if nargin >= 10 % DRM_DABNet(...,u_idx,y_idx)
                obj.u_idx = u_idx; obj.y_idx = y_idx;
                if nargin >= 12 % DRM_DABNet(...,u_idx,y_idx,u_name,y_name)
                    obj.u_name = u_name(u_idx); obj.y_name = y_name(y_idx);
                    if nargin >= 14 % DRM_DABNet(...,u_idx,y_idx,u_name,y_name,u_unit,y_unit)
                        obj.u_unit = u_unit(u_idx); obj.y_unit = y_unit(y_idx);
                        if nargin >= 15 % DRM_DABNet(...,u_idx,y_idx,u_name,y_name,u_unit,y_unit,t_unit)
                            obj.t_unit = t_unit;
                        end
                    end
                end
            end
            % instantiate x{i,j} and y{i}
            obj.x = cell(obj.ny,obj.nu); 
            % instantiate A_ann{i}, B_ann{i}
            obj.A_ann = cell(obj.ny,1); obj.B_ann = cell(obj.ny,1);
            % instantiate nx_ann{i}
            obj.nx_ann = zeros(obj.ny,1);
            % instantiate prediction derivative matrix, dX_P/ddU_M{i}, which is again a cell of {P,M}
            obj.dxannP_dduM = cell(obj.ny,1);
            
            % define A_ann, B_ann, nx_ann
            Ai = cell(obj.nu,obj.nu); Bi = cell(obj.nu,obj.nu);
            for iy = 1:obj.ny
                for i = 1:obj.nu
                    for j = 1:obj.nu
                        if i == j
                            Ai{i,j} = obj.A{iy,i};
                            Bi{i,j} = obj.B{iy,i}./obj.u_sigma(i);
                        else
                            Ai{i,j} = zeros(size(obj.A{iy,i},1),size(obj.A{iy,j},2));
                            Bi{i,j} = zeros(size(obj.B{iy,i},1),size(obj.B{iy,j},2));
                        end
                    end
                end
                A_ann_ = cell2mat(Ai); B_ann_ = cell2mat(Bi);
                obj.nx_ann(iy) = length(A_ann_);
                obj.A_ann{iy} = sparse(A_ann_); obj.B_ann{iy} = sparse(B_ann_);
                %obj.A_ann{iy} = cell2mat(Ai); obj.B_ann{iy} = cell2mat(Bi);
                %obj.nx_ann(iy) = length(obj.A_ann{iy});
            end
        end
        
        function setupStateEstimation(obj,pars)
        if nargin == 1
            fprintf('-- No State Estimation parameters were provided -- \n');
        else
            obj.fStateEstimation = true;
            if isfield(pars,'FilterMethod')
                obj.FilterMethod = pars.FilterMethod;
            end
            
            if isfield(pars,'SEFormulation')
                obj.SEFormulation = pars.SEFormulation;
            end
            
            if isfield(pars,'d_idx')
                obj.d_idx = pars.d_idx;
            else
                obj.d_idx = 1:obj.nu;
            end
            obj.nd = length(obj.d_idx);

            if isfield(pars,'qd') % vector qd is E(wd); Qd = diag(qd.^2)
                obj.qd = reshape(pars.qd,[],1);
                if length(obj.qd) == length(obj.d_idx)
                    obj.Qd = diag(obj.qd.^2);
                elseif length(obj.qd) == 1
                    obj.Qd = obj.qd^2 * eye(length(obj.d_idx));
                else
                    disp('Warning: Incorrect dimensions for Qd provided. Default SD value 0.1 will be used for all disturbance variables...');
                    obj.Qd = 0.01 * eye(length(obj.d_idx));
                end
            else
                obj.Qd = 0.01 * eye(length(obj.d_idx));
            end
            obj.Qd = sparse(obj.Qd);
            
            if isfield(pars,'IncludeStates') && isa(pars.IncludeStates,'logical')
                obj.fIncludeStates = pars.IncludeStates;
            end
            
            if isfield(pars,'qx') % scalar qx is E(wx_i), common for all x; Qx = diag([qx^2 qx^2 .... qx^2])
                obj.qx = reshape(pars.qx,[],1); obj.qx = obj.qx(1);
                obj.Qx = obj.qx^2 * eye(sum(obj.nx_ann));
            else
                obj.Qx = 0.01 * eye(sum(obj.nx_ann));
            end
            obj.Qx = sparse(obj.Qx);
            
            if isfield(pars,'r') % vector r is E(v); R = diag(r.^2)
                obj.r = reshape(pars.r,[],1);
                if length(obj.r) == obj.ny
                    obj.R = diag(obj.r.^2);
                elseif length(obj.r) == 1
                    obj.R = obj.r^2 * eye(obj.ny);
                else
                    disp('Warning: Incorrect dimensions for R provided. Default SD value 1 will be used for all measured outputs...');
                    obj.R = eye(obj.ny);
                end
            else
                obj.R = eye(obj.ny);
            end
            obj.R = sparse(obj.R);
            
            switch lower(obj.FilterMethod)
                case 'ekf'
                    F = cell(obj.ny+1,obj.ny+1);
                    for i = 1:obj.ny
                        for j = 1:obj.ny
                            if i == j
                                F{i,j} = obj.A_ann{i};
                            else
                                F{i,j} = zeros(size(obj.A_ann{i},1),size(obj.A_ann{j},2));
                            end
                        end
                    end
                    j = obj.ny + 1;
                    for i = 1:obj.ny
                        F{i,j} = obj.B_ann{i}(:,obj.d_idx);
                    end
                    i = obj.ny + 1;
                    for j = 1:obj.ny
                        F{i,j} = zeros(length(obj.d_idx),size(obj.A_ann{j},2));
                    end
                    i = obj.ny + 1; j = obj.ny + 1;
                    F{i,j} = eye(length(obj.d_idx));
                    obj.KF.F = sparse(cell2mat(F));
            end
        end
        end

        function initialize(obj,u)
            % Initialize the state and output variables based on input variables assuming steady-state condition
            u_scaled = (u - obj.u_mean)./obj.u_sigma; % scale steady-state input vector
            for iy = 1:obj.ny
                for j = 1:obj.nu
                    n = length(obj.A{iy,j});
                    % steady-state solution: x = Ax + Bu or x = inv(I-A)Bu
                    obj.x{iy,j} = (eye(n)-obj.A{iy,j})\obj.B{iy,j}*u_scaled(j);
                end
            end
            [obj.y,obj.dy_du] = obj.funcStateToOutput(obj.x);            
            % Initial State-Estimation Parameters
            if obj.fStateEstimation
                obj.d = zeros(obj.nd,1);
                x_aug = [obj.getANNStates(obj.x);obj.d];
                nx = length(x_aug) - obj.nd;
                if ~isempty(obj.d_idx)
                    switch lower(obj.SEFormulation)
                        case {'rsi','rside','random-step input','random-step input disturbance','rri','rride','random-ramp input','random-ramp input disturbance','p','periodic'}
                            if obj.fIncludeStates
                                obj.Q = blkdiag(obj.Qx, obj.Qd);
                            else
                                obj.Q = blkdiag(zeros(nx), obj.Qd);
                            end
%                         case {'rri','rride','random-ramp input','random-ramp input disturbance'}
%                             if obj.fIncludeStates
%                                 obj.Q = blkdiag(obj.Qx, zeros(obj.nd), obj.Qd);
%                             else
%                                 obj.Q = blkdiag(zeros(nx), zeros(obj.nd), obj.Qd);
%                             end
                        otherwise
                            disp('Warning: Incorrect SE Formulation specified. Using Random-Step Input Disturbance Estimation.');
                            if obj.fIncludeStates
                                obj.Q = blkdiag(obj.Qx, obj.Qd);
                            else
                                obj.Q = blkdiag(zeros(nx), obj.Qd);
                            end                            
                    end    
                else
                    obj.Q = obj.Qx; % Only State Disturbance
                end
                obj.Q = sparse(obj.Q);
                
                obj.Mk = x_aug; % initial mean state estimate
                obj.Pk = 0.001*eye(length(obj.Mk)); % initial state covariance estimate
                obj.SDk = abs(obj.r.*obj.y);
            end
        end

        function evalNextStep(obj,u,yp)
            if obj.fStateEstimation && nargin > 2
                switch lower(obj.FilterMethod)
                    case 'ekf'
                        % Prediction Step
                        obj.Mk = obj.ffcn(obj.Mk,u);
                        obj.Pk = obj.KF.F*obj.Pk*obj.KF.F' + obj.Q;
                        [obj.y,obj.KF.H] = obj.hfcn(obj.Mk);
                        % Correction Step
                        PxC_ = obj.Pk*obj.KF.H';
                        Sk = obj.KF.H*PxC_ + obj.R;
                        Lk = PxC_/Sk;
                        obj.Mk = obj.Mk + Lk*(yp-obj.y);
                        obj.Pk = (eye(length(obj.Mk))-Lk*obj.KF.H)*obj.Pk;
                        %obj.Mk = obj.Mk + obj.Lkf*(yp-obj.y); % Using Steady-State Kalman Gain
                        for iy = 1:length(obj.y)
                            obj.SDk(iy) = sqrt(Sk(iy,iy));
                        end
                        obj.x = obj.getDecoupledStates(obj.Mk(1:end-obj.nd));
                        obj.d = obj.Mk(end-obj.nd+1:end);
                        %obj.x = obj.evalNextState(obj.x,u);
                        [obj.y,obj.dy_du] = obj.funcStateToOutput(obj.x);
                    case 'ukf'
                        try
                            [obj.Mk,obj.Pk] = ukf_predict1(obj.Mk,obj.Pk,@(x_aug)obj.ffcn(x_aug,u),obj.Q);
                        catch
                            %disp('Warning: Issue with Cholesky Decomposition (prediction step). Results may be inaccurate. Adjust Q/R or use EKF method instead.');
                            obj.Pk = nearestSPD(obj.Pk);
                            [obj.Mk,obj.Pk] = ukf_predict1(obj.Mk,obj.Pk,@(x_aug)obj.ffcn(x_aug,u),obj.Q);
                        end
                        try
                            [obj.Mk,obj.Pk,~,~,Sk] = ukf_update1(obj.Mk,obj.Pk,yp,@(x_aug)obj.hfcn(x_aug),obj.R);
                        catch
                            %disp('Warning: Issue with Cholesky Decomposition (updation step). Results may be inaccurate. Adjust Q/R or use EKF method instead.');
                            obj.Pk = nearestSPD(obj.Pk);
                            %obj.Pk = obj.removeCrossTerms(obj.Pk);
                            [obj.Mk,obj.Pk,~,~,Sk] = ukf_update1(obj.Mk,obj.Pk,yp,@(x_aug)obj.hfcn(x_aug),obj.R);
                        end
                        for iy = 1:obj.ny
                            obj.SDk(iy) = sqrt(Sk(iy,iy));
                        end
                        obj.x = obj.getDecoupledStates(obj.Mk(1:end-obj.nd));
                        obj.d = obj.Mk(end-obj.nd+1:end);
                        %obj.x = obj.evalNextState(obj.x,u);
                        [obj.y,obj.dy_du] = obj.funcStateToOutput(obj.x);
                end
            else % Additive Output Disturbance Only
                obj.x = obj.evalNextState(obj.x,u);
                [obj.y,obj.dy_du] = obj.funcStateToOutput(obj.x);
            end
        end
        
        function x_aug = ffcn(obj,x_aug,u)
            x_ = x_aug(1:end-obj.nd);
            d_ = x_aug(end-obj.nd+1:end);
            if ~isempty(obj.d_idx)
                u(obj.d_idx) = u(obj.d_idx) + d_;
            end
            x_ = obj.evalNextState(x_,u);
            x_aug = [x_; d_];
        end
        
        function [y,H] = hfcn(obj,x_aug)
            x_ = x_aug(1:end-obj.nd);
            % u not needed, so no need to extract d
            if nargout >= 2 % needed in EKF (derivative of h-function)
                y = zeros(obj.ny,1);
                H = zeros(obj.ny,length(x_aug));
                for iy = 1:obj.ny
                    x_ann = obj.getDecoupledStates(x_,iy);
                    [y_scaled,dy_dxann] = obj.ANN(iy).predict(x_ann);
                    y(iy) = y_scaled*obj.y_sigma(iy) + obj.y_mean(iy);
                    i = sum(obj.nx_ann(1:iy-1));
                    n = obj.nx_ann(iy);
                    H(iy,i+1:i+n) = dy_dxann;
                end
                H = sparse(H);
            else
                y = obj.funcStateToOutput(x_);
            end
        end

        function [Y,dY_ddU,d2Y_ddU] = predictNextStep(obj,u,varargin)
            % Predict response (not update DRM) using defined input(s) (general)
            % Usage: [y_P,dy_du,d2y_ddu] = predictNextStep(u,P)
            %   where, size(u,2) ~ M
            P = 1;
            if nargin > 2
                P = varargin{1};
            end
            M = size(u,2);
            if  nargout >=2
                obj.generatePredictionMatrices(P,M);
            end
            if M < P
                u = [u repmat(u(:,end),1,P-M)];
            else
                u = u(:,1:P);
            end
            Y = repmat(obj.y,1,P);
            if nargout >= 2
                dY_ddU = cell(P,M);
                dY_ddU(:) = {zeros(obj.ny,obj.nu)};
            end
            if nargout >=3
                d2Y_ddU = cell(obj.ny,P); d2Y_ddU(:) = {cell(M,M)};
                for p = 1:P
                    for iy = 1:obj.ny
                        d2Y_ddU{iy,p}(:) = {zeros(obj.nu)};
                    end
                end
            end

            x_ = obj.x;
            for p = 1:P
                % predict states for p-th step
                u_ = u(:,p);
                if obj.fStateEstimation
                    if ~isempty(obj.d_idx)
                        u_(obj.d_idx) = u_(obj.d_idx) + obj.d;
                    end
                end
                x_ = obj.evalNextState(x_,u_);
                % evaluate outputs / derivatives at predicted states
                for iy = 1:obj.ny 
                    x_ann = obj.getANNStates(x_,iy);
                    if nargout >= 3 % Jacobian + Hessian Matrix
                        [ym_scaled,dy_dxann,d2y_dxann2] = obj.ANN(iy).predict(x_ann);
                    elseif nargout >= 2  % Jacobian Matrix only
                        [ym_scaled,dy_dxann] = obj.ANN(iy).predict(x_ann);
                    else
                        [ym_scaled] = obj.ANN(iy).predict(x_ann);
                    end
                    Y(iy,p) = ym_scaled*obj.y_sigma(iy) + obj.y_mean(iy); % de-scale output-vector
                    
                    if nargout >= 2
                        for m = 1:min(p,M)
                            dY_ddU{p,m}(iy,:) = (dy_dxann.*obj.y_sigma(iy)) * obj.dxannP_dduM{iy}{p,m};
                        end
                    end
                    
                    if nargout >=3
                        for m1 = 1:M
                            for m2 = 1:m1
                                d2Y_ddU{iy,p}{m1,m2} = obj.dxannP_dduM{iy}{p,m1}.' * (d2y_dxann2.*obj.y_sigma(iy)) * obj.dxannP_dduM{iy}{p,m2};
                                if m1~=m2
                                    d2Y_ddU{iy,p}{m2,m1} = d2Y_ddU{iy,p}{m1,m2}.';
                                end
                            end
                        end
                    end
                    
                end % iy
            end % p
        end % function
        
        function generatePredictionMatrices(obj,P,M)
            if obj.derivativeP~=P || obj.derivativeM~=M
                %disp('Inside generatePredictionMatrices');
                %tstart = tic;
                J = cell(P,M);
                for iy = 1:obj.ny
                    %J(:) = {zeros(size(obj.B_ann{iy},1),size(obj.B_ann{iy},2))}; % takes care of the left triangular zeros                    
                    J(:) = {sparse(size(obj.B_ann{iy},1),size(obj.B_ann{iy},2))}; % takes care of the left triangular zeros
                    % fill right triangular elements, one diagonal at a time
                    %pow = eye(size(obj.A_ann{iy})); sum = zeros(size(obj.B_ann{iy}));
                    pow = speye(length(obj.A_ann{iy})); sum = sparse(size(obj.B_ann{iy},1),size(obj.B_ann{iy},2));
                    for k = 1:P
                        sum = sum + pow*obj.B_ann{iy};
                        for j = 1:min(P-k+1,M)
                            J{j+k-1,j} = sum;
                        end
                        pow = obj.A_ann{iy}*pow;
                    end
                    obj.dxannP_dduM{iy} = J;
                end
                obj.derivativeP = P; obj.derivativeM = M;
                %fprintf('"Generate Prediction Matrix" duration = %d sec\n',toc(tstart));
                %pause
            end
        end
        
        function x = evalNextState(obj,x,u)
            if ~iscell(x)
                x = obj.getDecoupledStates(x);
                fANNFormRequired = true;
            else
                fANNFormRequired = false;
            end
            u_scaled = (u(:,1) - obj.u_mean)./obj.u_sigma; % scale input vector
            for iy = 1:obj.ny
                for j = 1:obj.nu
                    x{iy,j} = obj.A{iy,j}*x{iy,j} + obj.B{iy,j}*u_scaled(j);
                end
            end
            if fANNFormRequired
                x = obj.getANNStates(x);
            end
        end
        
        function [y,dy_du] = funcStateToOutput(obj,x)
            y = zeros(obj.ny,1);
            if nargout >= 2
                dy_du = zeros(obj.ny,obj.nu);
            end
            for iy = 1:obj.ny
                if ~iscell(x)
                    x_ann = obj.getDecoupledStates(x,iy);
                else
                    x_ann = obj.getANNStates(x,iy);
                end
                if nargout >= 2
                    [y_scaled,dy_dxann] = obj.ANN(iy).predict(x_ann);
                else
                    y_scaled = obj.ANN(iy).predict(x_ann);
                end
                y(iy) = y_scaled*obj.y_sigma(iy) + obj.y_mean(iy);
                if nargout >= 2
                    dy_du(iy,:) = dy_dxann*obj.B_ann{iy};
                end
            end
        end
        
        function x = getDecoupledStates(obj,x_ann,iy)
            if nargin > 2 % decouple states corresponding to output iy into x_ann
                i = sum(obj.nx_ann(1:iy-1));
                n = obj.nx_ann(iy);
                x = x_ann(i+1:i+n);
            elseif nargin > 1 % decouple all states into x{iy,j}
                x = cell(obj.ny,obj.nu);
                i = 0;
                for iy = 1:obj.ny
                    for j = 1:obj.nu
                        n = length(obj.A{iy,j});
                        x{iy,j} = x_ann(i+1:i+n);
                        i = i + n;
                    end
                end
            end
        end
        
        function x_ann = getANNStates(obj,x,iy)
            if nargin > 2 % populate states corresponding to output iy
                x_ann = zeros(obj.nx_ann(iy),1);
                i = 0;
                for j = 1:obj.nu
                    n = length(obj.A{iy,j});
                    x_ann(i+1:i+n) = x{iy,j};
                    i = i + n;
                end
            elseif nargin > 1 % populate all states (corresponding to all outputs)
                x_ann = zeros(sum(obj.nx_ann),1);
                i = 0;
                for iy = 1:obj.ny
                    for j = 1:obj.nu
                        n = length(obj.A{iy,j});
                        x_ann(i+1:i+n) = x{iy,j};
                        i = i + n;
                    end
                end
            end
        end

        function A = removeCrossTerms(obj,A)
            I = eye(length(A));
            i = 0;
            for iy = 1:obj.ny
                for j = 1:obj.nu
                    n = length(obj.A{iy,j});
                    I(i+1:i+n,i+1:i+n) = ones(n);
                    i = i + n;
                end
            end
            A = A.*I;
        end
        
        function evalNextStateRelatedToOutputs(obj,u,iy)
        % calculate x{i,j} and then update x_ann{} related to defined outputs
        % iy is vector of indices of related output variables
        % this function is called by Kalman filter's predict function
            u_scaled = (u - obj.u_mean)./obj.u_sigma; % Scale input vector
            niy = size(iy,1);
            for i = 1:niy
                iyi = iy(i);
                for j = 1:obj.nu
                    obj.x{iyi,j} = obj.A{iyi,j}*obj.x{iyi,j} + obj.B{iyi,j}*u_scaled(j); % linear marching of states: x[k+1] = Ax[k] + Bu[k]
                end
                obj.populateAnnStatesFromDecoupledStates(iyi);
            end
        end

        function y_out = evalRelatedOutputsFromAnnStates(obj,iy)
        % calculate related outputs using the related states in x_ann{} without touching the decoupled state variables x{i,j}
        % iy is vector of indices of related output variables
            niy = size(iy,1);
            y_out = zeros(niy,1);
            for i = 1:niy
                [ym_scaled] = obj.predictSingleOutput(iy(i)); % nonlinear mapping of states to outputs through neural-network
                y_out(i) = ym_scaled*obj.y_sigma(iy(i)) + obj.y_mean(iy(i)); % de-scale output-vector
            end
        end

        function populateAnnStatesFromDecoupledStates(obj,i)
        % set ANN states x_ann{i} using decoupled states in x{i,:}
        % i is the index of output
            ifirst = 1;
            for j = 1:obj.nu
                ilast = ifirst + length(obj.A{i,j}) - 1;
                obj.x_ann{i}(ifirst:ilast) = obj.x{i,j}(1:end);
                ifirst = ilast + 1;
            end
        end

        function populateDecoupledStatesFromAnnStates(obj,i)
        % set decoupled states x{i,:} using ANN states x_ann{i}
        % i is the index of output
            ifirst = 1;
            for j = 1:obj.nu
                ilast = ifirst + length(obj.A{i,j}) - 1;
                obj.x{i,j}(1:end) = obj.x_ann{i}(ifirst:ilast);
                ifirst = ilast + 1;
            end
        end

        function x_all = getAnnStatesRelatedToOutputs(obj,iy)
        % get states related to outputs in x_ann{}
        % x_all is a vector of all related state variables
        % iy is an array of indices of related outputs
            niy = size(iy,1);
            nstate = 0;
            for i = 1:niy
               nstate = nstate + obj.nx_ann(iy(i));
            end
            x_all = zeros(nstate,1);
            ifirst = 1;
            for i = 1:niy
                ilast = ifirst + obj.nx_ann(iy(i)) - 1;
                x_all(ifirst:ilast) = obj.x_ann{iy(i)}(1:end);
                ifirst = ilast + 1;
            end
        end

        function setAnnStatesRelatedToOutputs(obj,x_all,iy)
        % set states related to outputs in x_ann{}
        % x_all is a vector of all related state variables
        % iy is an array of indices of related outputs
            niy = size(iy,1);
            ifirst = 1;
            for i = 1:niy
                ilast = ifirst + obj.nx_ann(iy(i)) - 1;
                obj.x_ann{iy(i)}(1:end) = x_all(ifirst:ilast);
                ifirst = ilast + 1;
            end
        end

        function setAnnAndDecoupledStatesRelatedToOutputs(obj,x_all,iy)
        % set states related to ouputs in x_ann{} and x{i,j}
        % x_all is a vector of all related state variables
        % iy is an array of indices of related outputs
            niy = size(iy,1);
            ifirst = 1;
            for i = 1:niy
                ilast = ifirst + obj.nx_ann(iy(i)) - 1;
                obj.x_ann{iy(i)}(1:end) = x_all(ifirst:ilast);
                obj.populateDecoupledStatesFromAnnStates(iy(i));
                ifirst = ilast + 1;
            end
        end

        function [y,dy_dxann] = predictSingleOutput(obj,i)
        % do neural network mapping from state variables to scaled output
        % i is the output index
            if nargout >= 2
                [y,dy_dxann] = obj.ANN(i).predict(obj.x_ann{i});
            else
                [y] = obj.ANN(i).predict(obj.x_ann{i});
            end
        end

        function nx = getNumberOfStatesRelatedToOutputs(obj,iy)
        % get number of states related to specified outputs
        % iy is vector of indices of related output variables
            niy = size(iy,1);
            nx = 0;
            for i = 1:niy
                nx = nx + obj.nx_ann(iy(i)); 
            end
        end

        function icross = getCrossIndices(obj,iy)
        % returned icross is an augmented matrix with 1 for sub-matrices of x_ann{} and 0 in anywhere else, i.e. 0 for decoupled outputs
        % iy is vector of indices of related output variables
            niy = size(iy,1);
            nx = obj.getNumberOfStatesRelatedToOutputs(iy);
            icross = zeros(nx,nx);
            ifirst = 1;
            for i = 1:niy
                ilast = ifirst + obj.nx_ann(iy(i)) - 1;
                for ii = ifirst:ilast
                    for jj = ifirst:ilast
                        icross(ii,jj) = 1;
                    end
                end
                ifirst = ilast + 1;
            end
        end

        function icross = getCrossIndices2(obj,iy)
        % returned icross is an augmented matrix with 1 for sub-matrices of A{i,j} and 0 in anywhere else
        % iy is vector of indices of related output variables
            niy = size(iy,1);
            nx = obj.getNumberOfStatesRelatedToOutputs(iy);
            icross = zeros(nx,nx);
            ifirst = 1;
            for i = 1:niy
                for j = 1:obj.nu
                    ilast = ifirst + length(obj.A{iy(i),j}) - 1;
                    for ii = ifirst:ilast
                        for jj = ifirst:ilast
                            icross(ii,jj) = 1;
                        end
                    end
                    ifirst = ilast + 1;
                end
            end
        end
        
    end % methods
    
end % class
