classdef Quad_Tank < handle

% Quadruple-Tank Process
%
% Johansson, K.H., "The Quadruple-Tank Process: A Multivariable Laboratory
%   Process with an Adjustable Zero," IEEE Trans. Cont. Sys. Tech., 8(3), 2000
%
% Code written and developed by:
% Priyadarshi Mahapatra, PhD
% URS Corporation / National Energy Technology Laboratory
% 01 August 2014
% 
% States x correspond to [h1; h2; h3; h4]
% Outputs y correspond to [h1; h2; h3; h4]
%   Note - The above ref. considers voltage measurement corresponding to
%   heights as outputs [kc*h1; kc*h2] - This is modified to increase the
%   number of outputs for state-estimation studies and for brevity.

    properties
        A1 = 28;        % Cross-Sectional Area of Tank-1 [cm^2]
        A2 = 32;        % Cross-Sectional Area of Tank-2 [cm^2]
        A3 = 28;        % Cross-Sectional Area of Tank-3 [cm^2]
        A4 = 32;        % Cross-Sectional Area of Tank-4 [cm^2]
        
        a1 = 0.071;     % Orifice Area of Tank-1 [cm^2]
        a2 = 0.057;     % Orifice Area of Tank-2 [cm^2]
        a3 = 0.071;     % Orifice Area of Tank-3 [cm^2]
        a4 = 0.057;     % Orifice Area of Tank-4 [cm^2]
        
        %kc = 0.5;       % Elec. Coefficient for Height Measurement [V/cm]
        
        g = 981;        % Acceleration of Gravity [cm/s^2]
        
        k1;             % Valve Elec. Coefficient [cm^3/V/s]
        k2;             % Valve Elec. Coefficient [cm^3/V/s]
        
        gamma1;         % Split Fraction - T1/(T1+T4)
        gamma2;         % Split Fraction - T2/(T2+T3)
        
        v1;             % Voltage to Pump-1 [V]
        v2;             % Voltage to Pump-2 [V]
        
        h1;             % Height of Fluid in Tank-1 [cm]
        h2;             % Height of Fluid in Tank-2 [cm]
        h3;             % Height of Fluid in Tank-3 [cm]
        h4;             % Height of Fluid in Tank-4 [cm]

        NMP;            % Denotes Non-Minimum Phase State
        
        u_name = {'v_1','v_2'};
        y_name = {'h_1','h_2','h_3','h_4'};
        u_unit = {'Volt','Volt'};
        y_unit = {'cm','cm','cm','cm'};
        t_unit = 'sec';
        
        x; y; u;
        
        x_nominal; y_nominal; u_nominal;
        
        noise_u = 0.0005*1; % Input noise SD
        noise_y = 0.0005*1; % Output noise SD
    end
    
    methods
        function initialize(obj,NMP)
            obj.NMP = NMP;
            
            if NMP == true
                % NON-MINIMUM PHASE (large inverse response, difficult to control)
                obj.v1 = 3.15; obj.v2 = 3.15;
                obj.k1 = 3.14;  obj.k2 = 3.29;
                obj.gamma1 = 0.43;  obj.gamma2 = 0.34;
                % Initial condition based on reference [Johansson, 2000]
%                 obj.h1 = 12.6; obj.h2 = 13.0;
%                 obj.h3 = 4.8;  obj.h4 = 4.9;
                % Initial condition based on exact steady-state solution
                x_guess = [12.6; 13.0; 4.8; 4.9]; % Initial guess for solving SS
                options = optimset('Display','off');
                obj.x = fsolve(@(x) obj.qt_ODE(x),x_guess,options);
                obj.h1 = obj.x(1); obj.h2 = obj.x(2);
                obj.h3 = obj.x(3); obj.h4 = obj.x(4);
            else
                % MINIMUM PHASE
                obj.v1 = 3; obj.v2 = 3;
                obj.k1 = 3.33;  obj.k2 = 3.35;
                obj.gamma1 = 0.7;  obj.gamma2 = 0.6;
                % Initial condition based on reference [Johansson, 2000]
%                 obj.h1 = 12.4; obj.h2 = 12.7;
%                 obj.h3 = 1.8;  obj.h4 = 1.4;
                % Initial condition based on exact steady-state solution
                x_guess = [12.4; 12.7; 1.8; 1.4]; % Initial guess for solving SS
                options = optimset('Display','off');
                obj.x = fsolve(@(x) obj.qt_ODE(x),x_guess,options);
                obj.h1 = obj.x(1); obj.h2 = obj.x(2);
                obj.h3 = obj.x(3); obj.h4 = obj.x(4);
            end
            
            obj.u = [obj.v1; obj.v2];
            obj.x = [obj.h1; obj.h2; obj.h3; obj.h4];
            %obj.y = obj.kc .* [obj.h1; obj.h2];
            obj.y = [obj.h1; obj.h2; obj.h3; obj.h4];
            
            obj.u_nominal = obj.u;
            obj.x_nominal = obj.x;
            obj.y_nominal = obj.y;
        end
        
        function evalNextStep(obj,u,dt)
            u = u + obj.noise_u.*obj.u_nominal.*randn(size(obj.u_nominal)); % <----- force some input noise
            obj.u = u;
            obj.v1 = u(1); obj.v2 = u(2);
            % Plant Simulation ('dummy' represents continuous data b/w t_k and t_k+n)
            [~,xdummy] = ode45(@(t,x)obj.qt_ODE(x),[0 dt],obj.x);
            obj.x = xdummy(end,:)'; % Extract last point for sampling
            obj.h1 = obj.x(1); obj.h2 = obj.x(2);
            obj.h3 = obj.x(3); obj.h4 = obj.x(4);
            %obj.y = obj.kc .* [obj.h1; obj.h2];
            obj.y = [obj.h1; obj.h2; obj.h3; obj.h4];
            obj.y = obj.y + obj.noise_y.*obj.y_nominal.*randn(size(obj.y_nominal)); % <----- force some measurement noise
        end
        
        function xdot = qt_ODE(obj,x)
            H1 = x(1); H2 = x(2);
            H3 = x(3); H4 = x(4);
            % Differential equations
            dH1dt = (-obj.a1*sqrt(2*obj.g*H1) + obj.a3*sqrt(2*obj.g*H3) + obj.gamma1*obj.k1*obj.v1) / obj.A1;
            dH2dt = (-obj.a2*sqrt(2*obj.g*H2) + obj.a4*sqrt(2*obj.g*H4) + obj.gamma2*obj.k2*obj.v2) / obj.A2;
            dH3dt = (-obj.a3*sqrt(2*obj.g*H3) + (1-obj.gamma2)*obj.k2*obj.v2) / obj.A3;
            dH4dt = (-obj.a4*sqrt(2*obj.g*H4) + (1-obj.gamma1)*obj.k1*obj.v1) / obj.A4;
            xdot = [dH1dt; dH2dt; dH3dt; dH4dt];
        end
    end
end    

