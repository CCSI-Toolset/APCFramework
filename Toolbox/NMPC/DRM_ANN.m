classdef DRM_ANN < handle
    
% Class-definition for artificial neural network used by D-RM Builder
%
% Written and developed by:
% Priyadarshi Mahapatra, URS Corporation / National Energy Technology Laboratory, 2013
% Jinliang Ma, URS Corporation / National Energy Technology Laboratory, 2013
%
% Last Modified: September 2013

    properties
        nx;                 %number of input variables
        ny;                 %number of output variables
        nlayer_hid;         %number of hidden layers, always 1
        nneuron_hid;        %number of neurons in hidden layer
        steepness_hid;      %steepness in hidden layers
        steepness_out;      %steepness in output layer
        iactivation_hid;    %activation function type in hidden layers
        iactivation_out;    %activation function type in output layer
        weight;             %weight vector
        mean_in;            %mean vector of input variables
        sigma_in;           %standard deviation vector of input variables
        mean_out;           %mean vector of output variables
        sigma_out;          %standard deviation vector of output variables
        func_count = 0;     %no of times predict function is called
    end
    
    methods
        function obj = DRM_ANN(ann)
            if nargin == 1
                obj.nx = ann.nx;
                obj.ny = ann.ny;
                obj.nlayer_hid = ann.nhid;
                obj.nneuron_hid = ann.nneuron_hid;
                obj.steepness_hid = ann.steepness_hidden;
                obj.steepness_out = ann.steepness_output;
                obj.iactivation_hid = ann.iactivation_hidden;
                obj.iactivation_out = ann.iactivation_output;
                obj.weight = ann.weight;
                obj.mean_in = ann.mean_in;
                obj.mean_out = ann.mean_out;
                obj.sigma_in = ann.sigma_in;
                obj.sigma_out = ann.sigma_out;
            end
        end
        
        function [y,dy_dx,d2y_dx2] = predict(obj,x)
            obj.func_count = obj.func_count + 1;
            Narg = nargout; 
            
            % Memory Allotment for Values, Jacobians and Hessians
            y = zeros(1,obj.ny);            
            if Narg >= 2
                if obj.ny > 1
                    dy_dx = cell(1,obj.ny);
                    dy_dx(:) = {zeros(1,obj.nx)};
                else
                    dy_dx = zeros(1,obj.nx);
                end
                if Narg >= 3
                    if obj.ny > 1
                        d2y_dx2 = cell(1,obj.ny);
                        d2y_dx2(:) = {zeros(obj.nx)};
                    else
                        d2y_dx2 = zeros(obj.nx);
                    end
                end
            end            
            
%             Narg = 2; % To force zero hessian [accurate Hessian is sometimes not needed]

            y1 = zeros(1,obj.nx);
            y2 = zeros(1,obj.nneuron_hid);
            y3 = zeros(1,obj.ny);
            if Narg >= 2
                dydt2 = zeros(1,obj.nneuron_hid);
                dydt3 = zeros(1,obj.ny);
                if Narg >=3
                    dydtdy2 = zeros(1,obj.nneuron_hid);
                    dydtdy3 = zeros(1,obj.ny);
                end                    
            end
            % normalize input x and assign to input layer neuron output
            for i=1:obj.nx
                y1(i) = (x(i)-obj.mean_in(i))/obj.sigma_in(i);
            end
            % assign weight index - 1
            iw_start2 = 1:obj.nneuron_hid;
            iw_start3 = 1:obj.ny;
            iconn = 0;
            for i=1:obj.nneuron_hid
                iw_start2(i) = iconn;
                iconn = iconn + obj.nx + 1;
            end
            for i=1:obj.ny
                iw_start3(i) = iconn;
                iconn = iconn + obj.nneuron_hid + 1;
            end
            % calculate hidden layer data
            for i=1:obj.nneuron_hid
                iw_start = iw_start2(i);
                sum = 0;
                for j=1:obj.nx
                    sum = sum + y1(j)*obj.weight(iw_start+j);
                end
                % bias neuron contribution
                sum = sum + obj.weight(iw_start+obj.nx+1);
                sum = sum*obj.steepness_hid;
                switch obj.iactivation_hid
                case 0
                    y2(i) = sum;
                    if Narg >= 2
                        dydt2(i) = obj.steepness_hid;
                        if Narg >= 3
                            dydtdy2(i) = 0;
                        end
                    end
                case 1
                    y2(i) = 1/(1 + exp(-2*sum));
                    if Narg >= 2
                        dydt2(i) = 2*obj.steepness_hid*y2(i)*(1-y2(i));
                        if Narg >= 3
                            dydtdy2(i) = 2*obj.steepness_hid*(1-2*y2(i));
                        end
                    end
                case 2
                    y2(i) = 2/(1+exp(-2*sum))-1;
                    if Narg >= 2
                        dydt2(i) = obj.steepness_hid*(1-y2(i)*y2(i));
                        if Narg >= 3
                            dydtdy2(i) = -2*obj.steepness_hid*y2(i);
                        end
                    end
                otherwise
                    y2(i) = sum;
                    if Narg >= 2
                        dydt2(i) = obj.steepness_hid;
                        if Narg >= 3
                            dydtdy2(i) = 0;
                        end
                    end
                end
            end
            % calculate output layer data
            for i=1:obj.ny
                iw_start = iw_start3(i);
                sum = 0;
                for j=1:obj.nneuron_hid
                    sum = sum + y2(j)*obj.weight(iw_start+j);
                end
                % bias neuron contribution
                sum = sum + obj.weight(iw_start+obj.nneuron_hid+1);
                sum = sum*obj.steepness_out;
                switch obj.iactivation_out
                case 0
                    y3(i) = sum;
                    if Narg >= 2
                        dydt3(i) = obj.steepness_out;
                        if Narg >= 3
                            dydtdy3(i) = 0;
                        end
                    end
                case 1
                    y3(i) = 1/(1 + exp(-2*sum));
                    if Narg >= 2
                        dydt3(i) = 2*obj.steepness_out*y3(i)*(1-y3(i));
                        if Narg >= 3
                            dydtdy3(i) = 2*obj.steepness_out*(1-2*y3(i));
                        end
                    end
                case 2
                    y3(i) = 2/(1+exp(-2*sum))-1;
                    if Narg >= 2
                        dydt3(i) = obj.steepness_out*(1-y3(i)*y3(i));
                        if Narg >= 3
                            dydtdy3(i) = -2*obj.steepness_out*y3(i);
                        end
                    end
                otherwise
                    y3(i) = sum;
                    if Narg >= 2
                        dydt3(i) = obj.steepness_out;
                        if Narg >= 3
                            dydtdy3(i) = 0;
                        end
                    end
                end
                % ------------------------- Values
                y(i) = y3(i)*obj.sigma_out(i) + obj.mean_out(i);
                % ------------------------- Jacobians
                if Narg >= 2
                    % calculate Jacobian for scaled input and output variables
                    for j=1:obj.nx
                        sum = 0;
                        for k=1:obj.nneuron_hid
                            sum = sum + obj.weight(iw_start3(i)+k)*obj.weight(iw_start2(k)+j)*dydt2(k);
                        end
                        if obj.ny > 1
                            dy_dx{i}(j) = sum*dydt3(i)*obj.sigma_out(i)/obj.sigma_in(j);
                        else
                            dy_dx(j) = sum*dydt3(i)*obj.sigma_out(i)/obj.sigma_in(j);
                        end                        
                    end % j
                % ------------------------- Hessians
                    if Narg >= 3
                        for j=1:obj.nx
                            for l=j:obj.nx
                                if obj.ny > 1
                                    term1 = dy_dx{i}(j)/dydt3(i)*dydtdy3(i)*dy_dx{i}(l);
                                else
                                    term1 = dy_dx(j)/dydt3(i)*dydtdy3(i)*dy_dx(l);
                                end
                                sum = 0;
                                for k=1:obj.nneuron_hid
                                    sum = sum + obj.weight(iw_start3(i)+k)*obj.weight(iw_start2(k)+j)*dydtdy2(k)*obj.weight(iw_start2(k)+l)*dydt2(k);
                                end
                                term2 = dydt3(i)*sum;
                                if obj.ny > 1
                                    d2y_dx2{i}(j,l) = (term1 + term2)*obj.sigma_out(i)/obj.sigma_in(j)/obj.sigma_in(l);
                                    if j~=l
                                        d2y_dx2{i}(l,j) = d2y_dx2{i}(j,l);
                                    end
                                else
                                    d2y_dx2(j,l) = (term1 + term2)*obj.sigma_out(i)/obj.sigma_in(j)/obj.sigma_in(l);
                                    if j~=l
                                        d2y_dx2(l,j) = d2y_dx2(j,l);
                                    end
                                end
                            end % l
                        end % j
                        
%                         % gradNumerical = jacobianest(@(x)obj.predict(x),x);
%                         options = optimset('Display','off','GradObj','on','MaxFunEvals',0);
%                         [~,valNumerical,~,~,gradNumerical,hessNumerical] = fminunc(@(x)obj.predict(x),x,options);
%                         hessNumerical = full(hessNumerical);
%                         fprintf('Norm of Hessian difference (ana vs. num): %d \n',norm(hessNumerical-d2y_dx2));
%                         keyboard; pause;
                        
                    end % Narg >= 3
                end % Narg >= 2
                
            end % i
            
            
        end % function
    end % methods
    
end % classdef
