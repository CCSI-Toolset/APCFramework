classdef DRM_SS < handle

% Class-definition for STATE-SPACE based dynamic reduced-order models
%
% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
%
% Revision List
% 2013/10 - Initial version containing constructor, initialize and
%           evalNextStep and predictNextStep functions.
% 2014/08 - Modified evalNextStep and predictNextStep functions to accept M control moves
% 2015/03 - Code revision for computational speed including
%           a. Avoid fixed derivatives (dyP_dduM) evaluations everytime
%              predictNextStep function is invoked except when P, M changes
%              by introducing derivativeP & derivativeM properties
%           b. Preallocation of memory

    properties
        A; B; C; D; K;
        x;              % State(s) [deviation]
        y;              % Outputs(s) [physical]
        dy_du;          % Jacobian Matrix
        Ts;             % Sampling Time
        uss; yss;       % Steady-state value(s)
        u_idx; y_idx;   % Input/Output indices from among plant's I/O ports
        u_name;         % Input variable's names/tags
        y_name;         % Output variable's names/tags
        u_unit;         % Input variable' units
        y_unit;         % Output variable' units
        t_unit;         % Time unit
        nu; ny; nx;
        % Variables related to state-estimation
        fStateEstimation = false;
        fIncludeStates = false;
        d_idx;          % indices of estimated disturbances among model inputs, integer array;
        nd;
        d;              % estimated disturbances
        Mk; Pk; SDk;
        KF;
        %Lkf;            % Kalman Filter Gain
        FilterMethod = 'KF'; % Choose between 'KF' or 'UKF'
        SEFormulation = 'Random-Step Input Disturbance';
        qd; qx; r;
        Qd; Qx;
        Q; R;
    end        
    properties (SetAccess = private)
        derivativeP = 0; derivativeM = 0;
        Sf;
    end
    
    methods
        function obj = DRM_SS(varargin)
            switch nargin
                case 1 % DRM_SS(IdentMod)
					IdentMod = varargin{1};
                    obj.A = IdentMod.A; obj.B = IdentMod.B;
                    obj.C = IdentMod.C; obj.D = IdentMod.D;
                    obj.K = IdentMod.K;
                    obj.x = IdentMod.x0; obj.Ts = IdentMod.Ts;
                case 5 % DRM_SS(A,B,C,D,Ts)
                    obj.A = varargin{1}; obj.B = varargin{2}; obj.C = varargin{3}; obj.D = varargin{4};
                    obj.Ts = varargin{5};					
                case 6 % DRM_SS(A,B,C,D,x0,Ts)
                    obj.A = varargin{1}; obj.B = varargin{2}; obj.C = varargin{3}; obj.D = varargin{4};
                    obj.x = varargin{5}; obj.Ts = varargin{6};
                case 7 % DRM_SS(A,B,C,D,x0,Ts,K)
                    obj.A = varargin{1}; obj.B = varargin{2}; obj.C = varargin{3}; obj.D = varargin{4};
                    obj.x = varargin{5}; obj.Ts = varargin{6}; obj.K = varargin{7};
                otherwise
                    disp('Incorrect state-space arguments provided. See documentation for correct usage.');
                    return;
            end
            obj.nu = size(obj.B,2); obj.ny = size(obj.C,1);
            obj.nx = size(obj.A,1); 
        end
        
%         function setupStateEstimation(obj,FilterMethod,q,r,d_idx)
%             obj.fStateEstimation = true;
%             obj.FilterMethod = FilterMethod;
%             obj.q_spec = q; obj.r_spec = r;
%             if nargin > 4
%                 obj.d_idx = d_idx;
%             else
%                 obj.d_idx = 1:length(obj.u_idx);
%             end
%             obj.nd = length(obj.d_idx);            
%         end

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
            
            if isfield(pars,'IncludeStates') && isa(pars.IncludeStates,'logical')
                obj.fIncludeStates = pars.IncludeStates;
            end
            
            if isfield(pars,'d_idx')
                obj.d_idx = pars.d_idx;
            else
                obj.d_idx = 1:obj.nu;
            end
            obj.nd = length(obj.d_idx);

            if isfield(pars,'qd') % vector qd is E(wd); Qd = diag(qd.^2)
                obj.qd = reshape(pars.qd,[],1);
            else
                obj.qd = 0.1;
            end
                        
            if isfield(pars,'qx') % scalar qx is E(wx_i), common for all x; Qx = diag([qx^2 qx^2 .... qx^2])
                obj.qx = reshape(pars.qx,[],1);
                obj.qx = obj.qx(1);
            else
                obj.qx = 0.1;
            end
            
            if isfield(pars,'r') % vector r is E(v); R = diag(r.^2)
                obj.r = reshape(pars.r,[],1);
            else
                obj.r = 1;
            end
        end
        end

        function initialize(obj,u)
            obj.x = (eye(length(obj.A))-obj.A)\obj.B*(u(:,1) - obj.uss);
            obj.y = obj.funcStateToOutput(obj.x,u);
            obj.dy_du = obj.C*obj.B + obj.D;
            if obj.fStateEstimation
                obj.d = zeros(obj.nd,1);              
                x_aug = [obj.x;obj.d];
                % Disturbance Covariance Matrix
                if length(obj.qd) == length(obj.d_idx)
                    obj.Qd = diag(obj.qd.^2);
                elseif length(obj.qd) == 1
                    obj.Qd = obj.qd^2 * eye(length(obj.d_idx));
                else
                    disp('Warning: Incorrect dimensions for Qd provided. Default SD value 0.1 will be used for all disturbance variables...');
                    obj.Qd = 0.1^2 * eye(length(obj.d_idx));
                end
                % State Covariance Matrix
                obj.Qx = obj.qx^2 * eye(obj.nx);
                % Measurement Covariance Matrix
                if length(obj.r) == obj.ny
                    obj.R = diag(obj.r.^2);
                elseif length(obj.r) == 1
                    obj.R = obj.r^2 * eye(obj.ny);
                else
                    disp('Warning: Incorrect dimensions for R provided. Default SD value 1 will be used for all measured outputs...');
                    obj.R = eye(obj.ny);
                end
                % obj.R = sparse(obj.R);
                                
                if ~isempty(obj.d_idx)
                    switch lower(obj.SEFormulation)
                        case {'rsi','rside','random-step input','random-step input disturbance','rri','rride','random-ramp input','random-ramp input disturbance','p','periodic'}
                            if obj.fIncludeStates
                                obj.Q = blkdiag(obj.Qx, obj.Qd);
                            else
                                obj.Q = blkdiag(zeros(obj.nx), obj.Qd);
                            end
%                         case {'rri','rride','random-ramp input','random-ramp input disturbance'}
%                             if obj.fIncludeStates
%                                 obj.Q = blkdiag(obj.Qx, zeros(obj.nd), obj.Qd);
%                             else
%                                 obj.Q = blkdiag(zeros(obj.nx), zeros(obj.nd), obj.Qd);
%                             end
                        otherwise
                            disp('Warning: Incorrect SE Formulation specified. Using Random-Step Input Disturbance Estimation.');
                            if obj.fIncludeStates
                                obj.Q = blkdiag(obj.Qx, obj.Qd);
                            else
                                obj.Q = blkdiag(zeros(obj.nx), obj.Qd);
                            end                                                        
                    end
                else
                    obj.Q = obj.Qx; % Only State Disturbance
                end
                % obj.Q = sparse(obj.Q);

                
                obj.Mk = x_aug; % initial mean state estimate
                obj.Pk = 0.001*eye(length(obj.Mk)); % initial state covariance estimate
                obj.SDk = abs(obj.r.*obj.y);
                switch lower(obj.FilterMethod)
                    case 'kf'
                        Bd = obj.B(:,obj.d_idx);
                        Dd = obj.D(:,obj.d_idx);
                        obj.KF.F = [obj.A Bd; zeros(obj.nd,obj.nx) eye(obj.nd)];
                        %obj.KF.B_aug = [obj.B; zeros(obj.nd,nu)];
                        obj.KF.H = [obj.C Dd];
                        %obj.KF.D_aug = obj.D;
                        
                        %Bw = [zeros(obj.nx,obj.nd); eye(obj.nd)];
                        %Dw = zeros(obj.ny,obj.nd);
                        %Qkf = obj.q_spec^2*eye(obj.nd); Rkf = obj.r_spec^2*eye(obj.ny);
                        %sys = ss(obj.KF.F,[obj.KF.B_aug Bw],obj.KF.H,[obj.KF.D_aug Dw],obj.Ts);
                        %[~,~,~,obj.Lkf] = kalman(sys,Qkf,Rkf);
                end
            end
        end

        function evalNextStep(obj,u,yp)
            if obj.fStateEstimation && nargin > 2
                switch lower(obj.FilterMethod)
                    case {'kf','ekf'}
                        % Prediction Step
                        obj.Mk = obj.ffcn(obj.Mk,u);
                        obj.Pk = obj.KF.F*obj.Pk*obj.KF.F' + obj.Q;
                        obj.y = obj.hfcn(obj.Mk,u);
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
                        obj.x = obj.Mk(1:end-obj.nd);
                        obj.d = obj.Mk(end-obj.nd+1:end);
                        obj.y = obj.hfcn(obj.Mk,u);
                        
                    case 'ukf' % Note: same response as linear 'kf' - may delete this along with ffcn and hfcn
                        [obj.Mk,obj.Pk] = ukf_predict1(obj.Mk,obj.Pk,@(x_aug)obj.ffcn(x_aug,u),obj.Q);
                        [obj.Mk,obj.Pk,~,~,Sk] = ukf_update1(obj.Mk,obj.Pk,yp,@(x_aug)obj.hfcn(x_aug,u),obj.R);
                        %obj.Pk = obj.removeCrossTerms(obj.Pk);
                        %obj.Pk = nearestSPD(obj.Pk);
                        for iy = 1:length(obj.y)
                            obj.SDk(iy) = sqrt(Sk(iy,iy));
                        end
                        obj.x = obj.Mk(1:end-obj.nd);
                        obj.d = obj.Mk(end-obj.nd+1:end);
                        obj.y = obj.hfcn(obj.Mk,u);
                end
            else % No State Estimation
                obj.x = obj.evalNextState(obj.x,u);
                obj.y = obj.funcStateToOutput(obj.x,u);
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
        
        function y = hfcn(obj,x_aug,u)
            x_ = x_aug(1:end-obj.nd);
            d_ = x_aug(end-obj.nd+1:end);
            if ~isempty(obj.d_idx)
                u(obj.d_idx) = u(obj.d_idx) + d_;
            end
            y = obj.funcStateToOutput(x_,u);
        end
        
        function [Y,dY_ddU] = predictNextStep(obj,u,varargin)
            % Predict response (not update DRM) using defined input(s) (general)
            % Usage: [y_P,dy_du] = predictNextStep(u,P)
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
            x_ = obj.x;
            for p = 1:P
                u_ = u(:,p);
                if obj.fStateEstimation
                    if ~isempty(obj.d_idx)
                        u_(obj.d_idx) = u_(obj.d_idx) + obj.d;
                    end
                end
                x_ = obj.evalNextState(x_,u_);
                Y(:,p) = obj.funcStateToOutput(x_,u_);
                %Y_N(:,k) = obj.CA{k}*obj.x + obj.S{k}*(u_N(:,k) - obj.uss) + obj.yss; % Faster
            end
            
            if nargout >= 2
                dY_ddU = obj.Sf;
            end
        end
        
        function dyP_dduM = generatePredictionMatrices(obj,P,M)
            if obj.derivativeP~=P || obj.derivativeM~=M
                %disp('Inside generatePredictionMatrices');
                %tstart = tic;
                obj.Sf = cell(P,M);
                obj.Sf(:) = {zeros(size(obj.C,1),size(obj.B,2))}; % takes care of the left triangular zeros
                % fill right triangular elements, one diagonal at a time
                pow = eye(size(obj.A)); sum = zeros(size(obj.B));
                for k = 1:P
                    sum = sum + pow*obj.B;
                    for j = 1:min(P-k+1,M)
                        obj.Sf{j+k-1,j} = obj.C*sum + obj.D;
                    end
                    pow = obj.A*pow;
                end
                obj.derivativeP = P; obj.derivativeM = M;
                dyP_dduM = obj.Sf;
                %fprintf('"Generate Prediction Matrix" duration = %d sec\n',toc(tstart));
                %pause
            end
        end
        
        function x = evalNextState(obj,x,u)
            x = obj.A*x + obj.B*(u(:,1) - obj.uss);
        end
        
        function y = funcStateToOutput(obj,x,u)
            y = obj.C*x + obj.D*(u(:,1) - obj.uss) + obj.yss;
        end
        
    end % methods
    
end % class
