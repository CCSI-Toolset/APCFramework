classdef pH_Reactor < handle
% Two-tank pH Neutralization Process
% ----------------------------------
%
% Component A (acid) is acetic-acid (HAC)
% Component B (base) is sodium hydroxide (NaOH)
%
% Inputs u correspond to [Q1; Q3; Q5; Q7]
% States x correspond to [h1; alpha1; beta1; h2; alpha2; beta2]
% Outputs y correspond to [pH1; pH2] 
%     where, pH1 = f(alpha1,beta1) and pH2 = f(alpha2,beta2)
%
% Usage notes -
% 1. Create a pH_Reactor object, for example:
%        pH_Neut_Plant = pH_Reactor;
% 2. Initialize the reactor object, for example:
%        pH_Neut_Plant.Initialize([2.46; 30.9; 3.0; 30.9]); % [Q1;Q3;Q5;Q7]
% 3. For obtaining outputs (y) at next sample-time use evalNextStep as:
%        pH_Neut_Plant.evalNextStep(u,Ts);
%        y = pH_Neut_Plant.y;
%
% For detailed description of the mathematical formulation and concepts, 
% please refer to:
% G.B. Sentoni and L.T. Biegler, State-Space Nonlinear Process Modeling: 
%     Identification and Universality, AIChE Journal, Vol. 44, No. 10, 
%     October 1998

% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
% 20 Nov 2013

    properties
        A1 = 10;		% [m^2]
        A2 = 10;		% [m^2]
        cv1 = 113.1;	% [m^(2.5)/hr]
        cv2 = 226.2;	% [m^(2.5)/hr]

        h1;             % Tank-1 Level [m]
        alpha1;     	% Conc. of A in Tank-1 [kmol/m^3]
        beta1;      	% Conc. of B in Tank-1 [kmol/m^3]	
        h2;     		% Tank-2 Level [m]
        alpha2;         % Conc. of A in Tank-2 [kmol/m^3]
        beta2;          % Conc. of B in Tank-2 [kmol/m^3]	

        Q1;     		% Acid Vol. Flowrate to Tank-1 [m^3/hr]
        Q2 = 2.4;       % Acid Vol. Flowrate to Tank-1 [m^3/hr]
        Q3;     		% Base Vol. Flowrate to Tank-1 [m^3/hr]

        Q5;             % Acid Vol. Flowrate to Tank-2 [m^3/hr]
        Q6 = 1.86;		% Acid Vol. Flowrate to Tank-2 [m^3/hr]
        Q7;             % Base Vol. Flowrate to Tank-2 [m^3/hr]

        Q4;             % Vol. Flowrate out from Tank-1
        Q8;             % Vol. Flowrate out from Tank-2

        Ca = 0.32;      % Acid (HAC) Feed Concentration [kmol/m^3]
        Cb = 0.05;      % Base (NaOH) Feed Concentration [kmol/m^3]
        
        pH1;            % pH in Tank-1
        pH2;            % pH in Tank-2

        u_name = {'Q_1'; 'Q_3'; 'Q_5'; 'Q_7'};
        y_name = {'pH_1'; 'pH_2'};
        u_unit = {'m^3/hr'; 'm^3/hr'; 'm^3/hr'; 'm^3/hr'};
        y_unit = {''; ''};
        t_unit = 'hr';
        
        x; y; u;
        
        x_nominal; y_nominal; u_nominal;
        
        noise_u = [0; 0; 0; 0]; % Input noise SD
        noise_y = [0; 0]; % Output noise SD
    end
    
    methods
        function initialize(obj,u)
            obj.u = u;
            obj.Q1 = u(1); obj.Q3 = u(2); obj.Q5 = u(3); obj.Q7 = u(4);
            
            % Define initial conditions for high-fidelity model or 'plant' @SS            
            x_guess = [0.1; 0.04; 0.04; 0.1; 0.04; 0.04]; % Initial guess for solving SS
            options = optimset('Display','off');
            obj.x = fsolve(@(x) obj.pH_ODE(x),x_guess,options);
            
            obj.h1 = obj.x(1); obj.alpha1 = obj.x(2); obj.beta1 = obj.x(3);
            obj.h2 = obj.x(4); obj.alpha2 = obj.x(5); obj.beta2 = obj.x(6);
            obj.Q4 = obj.cv1*sqrt(obj.h1); obj.Q8 = obj.cv2*sqrt(obj.h2);
            
            log_C_guess = [-7;-7;-3.5;-1.5];
            log_C = fsolve(@(x) obj.pH_AE(x,[obj.alpha1;obj.beta1]),log_C_guess,options);
            obj.pH1 = -log_C(1);
            log_C = fsolve(@(x) obj.pH_AE(x,[obj.alpha2;obj.beta2]),log_C_guess,options);
            obj.pH2 = -log_C(1);
            obj.y = [obj.pH1;obj.pH2]; % Initial plant-outputs
            
            obj.u_nominal = obj.u;
            obj.x_nominal = obj.x;
            obj.y_nominal = obj.y;
        end
        
        function evalNextStep(obj,u,dt)
            u = u + obj.noise_u.*obj.u_nominal.*randn(size(obj.u_nominal)); % <----- force some input noise
            obj.u = u;
            obj.Q1 = u(1); obj.Q3 = u(2); obj.Q5 = u(3); obj.Q7 = u(4);
            
            % Plant Simulation ('dummy' represents continuous data b/w t_k and t_k+n)
            [~,xdummy] = ode45(@(t,x) obj.pH_ODE(x),[0 dt],obj.x);
            obj.x = xdummy(end,:)'; % Extract last point for sampling
            
            obj.h1 = obj.x(1); obj.alpha1 = obj.x(2); obj.beta1 = obj.x(3);
            obj.h2 = obj.x(4); obj.alpha2 = obj.x(5); obj.beta2 = obj.x(6);
            obj.Q4 = obj.cv1*sqrt(obj.h1); obj.Q8 = obj.cv2*sqrt(obj.h2);
            
            options = optimset('Display','off');
            log_C_guess = [-7;-7;-3.5;-1.5];
            log_C = fsolve(@(x) obj.pH_AE(x,[obj.alpha1;obj.beta1]),log_C_guess,options);
            obj.pH1 = -log_C(1);
            log_C = fsolve(@(x) obj.pH_AE(x,[obj.alpha2;obj.beta2]),log_C_guess,options);
            obj.pH2 = -log_C(1);
            
            obj.y = [obj.pH1;obj.pH2];
            obj.y = obj.y + obj.noise_y.*obj.y_nominal.*randn(size(obj.y_nominal)); % <----- force some measurement noise
        end
        
        function xdot = pH_ODE(obj,x)
            h1 = x(1); alpha1 = x(2); beta1 = x(3);
            h2 = x(4); alpha2 = x(5); beta2 = x(6);
            Q1 = obj.Q1; Q2 = obj.Q2; Q3 = obj.Q3;
            Q5 = obj.Q5; Q6 = obj.Q6; Q7 = obj.Q7;
            Q4 = obj.cv1*sqrt(h1); Q8 = obj.cv2*sqrt(h2);
            Ca = obj.Ca; Cb = obj.Cb;

            % Differential equations
            dh1dt = (Q1 + Q2 + Q3 - Q4)/obj.A1;
            dalpha1dt = (Q1*Ca + Q2*Ca - Q4*alpha1)/h1/obj.A1;
            dbeta1dt = (Q3*Cb - Q4*beta1)/h1/obj.A1;

            dh2dt = (Q4 + Q5 + Q6 + Q7 - Q8)/obj.A2;
            dalpha2dt = (Q4*alpha1 + Q5*Ca + Q6*Ca - Q8*alpha2)/h2/obj.A2;
            dbeta2dt = (Q4*beta1 + Q7*Cb - Q8*beta2)/h2/obj.A2;

            xdot = [dh1dt; dalpha1dt; dbeta1dt; dh2dt; dalpha2dt; dbeta2dt];
        end
        
        function f = pH_AE(~,x,par)
            log_Ka = log10(1.8e-5);
            log_Kw = log10(1e-14);
            log_H = x(1); log_OH = x(2); log_HAC = x(3); log_AC = x(4);
            alpha = par(1); beta = par(2);

            f(1) = log_H + log_OH - log_Kw;
            f(2) = log_H + log_AC - log_HAC - log_Ka;
            f(3) = 10^log_HAC + 10^log_AC - alpha;
            f(4) = 10^log_OH + 10^log_AC - 10^log_H - beta;
        end
    end
end    

