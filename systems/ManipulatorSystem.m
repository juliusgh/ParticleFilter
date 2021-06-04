classdef ManipulatorSystem < System
    %MANIPULATORSYSTEM "Ebener Manipulator"
    %   Simple mechanical test system
    properties
        G = 9.81 % gravity constant
        L1 = 1 % rod 1 length
        L2 = 1 % rod 2 length
        S1 = 1/2 % rod 1 center
        S2 = 1/2 % rod 2 center
        M1 = 1 % rod 1 mass
        M2 = 1 % rod 2 mass
        I1 = 1/12 % rod 1 inertia torque
        I2 = 1/12 % rod 2 inertia torque
        T0 = 0.05 % sampling time
        MaxStepSimulation = 0.001 % maximum step for simulation
        MaxStepPrediction = 0.01 % maximum step for prediction
    end
    
    methods
        function obj = ManipulatorSystem()
            %MANIPULATORSYSTEM Construct an instance of this class
            %   Simple mechanical test system
            obj.dim = 4; % system dimensions
            obj.dimz = 4; % measurement dimensions
            obj.dimu = 2; % control input dimensions
            obj.x0 = [pi/2;pi;0;0]; % initial condition
            Vtech = HyperRect([-0.002 0.002],[-0.002 0.002],[-0.002 0.002],[-0.002 0.002]);
            obj.Vsim = HyperRect([-0.0 0.0],[-0.0 0.0],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.V = obj.Vsim + Vtech; % system noise support
%             obj.W = HyperRect([-0.1 0.1],[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]); % measurement noise support
            obj.Vsim = HyperRect([-0.0 0.0],[-0.0 0.0],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.W = HyperRect([-0.05 0.05],[-0.05 0.05],[-0.05 0.05],[-0.05 0.05]); % measurement noise support
            obj.Wsim = obj.W;
            %obj.Vsim = HyperRect([-0.0 0.0],[-0.0 0.0],[-0.0 0.0],[-0.0 0.0]); % system noise support
            obj.X0 = obj.x0 + obj.W; % initial noise support
            obj.uk = @(k) [0;0]; % fixed input function
        end
        
        function x_ = g(obj,x,u)
            %G process function with input
            %   x: current system state, u: control input
            if size(x,2) > 1
                x_ = NaN(size(x));
                for i = 1:size(x,2)
                    x_(:,i) = obj.g(x(:,i),u);
                end
            else
                options = odeset('MaxStep',obj.MaxStepSimulation);
                [~,y] = ode45(@(t,x) obj.f(t,x,u),[0 obj.T0],x,options);
                x_ = y(end,:)';
            end
        end
        
        function x_ = gpredict(obj,x,u)
            %G process function with input for prediction (less accurate)
            %   x: current system state, u: control input
            step = obj.MaxStepPrediction;
            n = ceil(obj.T0/step) - 1;
            for i=1:n
                x = obj.rk45(@(t,x) obj.f(t,x,u),0,step,x);
            end
            x_ = obj.rk45(@(t,x) obj.f(t,x,u),0,obj.T0-n*step,x);
        end
        
        function x_ = rk45(~,f,t,h,x)
			%RK45 Runge Kutta 4th order
			%   f: RHS, t: time, h: step, x: state
            k1 = f(t,x);
            k2 = f(t+h/2,x+h/2*k1);
            k3 = f(t+h/2,x+h/2*k2);
            k4 = f(t+h,x+h*k3);
            x_ = x + h/6*(k1+2*k2+2*k3+k4);
        end
        
        function x_ = f(obj,~,x,u)
			%F RHS of the ODE system
            [~,~,x3,x4] = unpack(x);
            x_ = [x3;x4;obj.n(x,u)];
        end
        
        function n = n(obj,x,u)
            %N non-linear process function part
            %   x: current system state
            [x1,x2,x3,x4] = unpack(x);
            [u1,u2] = unpack(u);
            s1 = sin(x1);
            s2 = sin(x2);
            s12 = sin(x1-x2);
            c12 = cos(x1-x2);
            a = obj.M1 * obj.S1^2 + obj.M2 * obj.L1^2 + obj.I1;
            b = obj.M2 .* obj.L1 .* obj.S2 .* c12;
            d = obj.M2 * obj.S2^2 + obj.I2;
            e = obj.M1 * obj.S1 + obj.M2 * obj.L1;
            n = 1./(a*d-b.^2) .* (obj.M2 .* obj.L1 .* obj.S2 .* s12 .* [-d.*(x4.^2)-b.*(x3.^2);b.*(x4.^2)+a.*(x3.^2)]...
                - obj.G .* [e.*d.*s1 - obj.M2.*obj.S2.*b.*s2;-e.*b.*s1 + obj.M2.*obj.S2.*a.*s2] + [d.*(u1-u2)-b.*u2;-b.*(u1-u2)+a.*u2]);
        end
        
        function z = h(obj,x)
            %H measurement function
            %   x: current system state
            %z = [eye(2) zeros(2)]*x;
            [x1,x2] = unpack(x);
            z = [obj.L1 .* sin(x1);...
                -obj.L1 .* cos(x1);...
                obj.L1 .* sin(x1) + obj.L2 .* sin(x2);...
                -obj.L1 .* cos(x1) - obj.L2 .* cos(x2)];
        end
        
        function u = uref(obj,xref)
            %UREF get desired control input for reference state
            %   x: reference system state, TODO: rewrite?
            xr = obj.f(0,xref,[0;0]);
            u = -xr(3:4,:);
        end
        
        function E = energy(obj,x)
            %ENERGY calculate total energy of the system, TODO: rewrite
            %   x: system states
            if nargin < 2
                x = obj.x;
            end
            Ekin = 0.5*(obj.M1*obj.S1^2 + obj.M2*obj.L1^2 + obj.I1).*x(3,:).^2 ... 
            +0.5*(obj.M2*obj.S2^2 + obj.I2).*x(4,:).^2 + obj.M2*obj.L1*obj.S2.*x(3,:).*x(4,:).*cos(x(1,:)-x(2,:));
            Epot = obj.G * ((obj.M1*obj.S1 + obj.M2*obj.L1) * (1 - cos(x(1,:))) + obj.M2*obj.S2*(1-cos(x(2,:))));
            E = Ekin + Epot;
        end
        
        function p = pos1(obj,x)
            if nargin < 2
                x = obj.x;
            end
            p = [obj.L1 * sin(x(1,:)); -obj.L1 * cos(x(1,:))];
        end
        
        function p = pos2(obj,x)
            if nargin < 2
                x = obj.x;
            end
            p = [obj.L1*sin(x(1,:))+obj.L2*sin(x(2,:));-obj.L1*cos(x(1,:))-obj.L2*cos(x(2,:))];
        end
        
        function p = getPos1(obj,t)
             k = floor(t / obj.T0) + 1;
             x = obj.x(:,k);
             p = obj.pos1(x);
        end
        
        function p = getPos2(obj,t)
             k = floor(t / obj.T0) + 1;
             x = obj.x(:,k);
             p = obj.pos2(x);
        end
        
        function g = plotPos(obj,t)
            p1 = obj.getPos1(t);
            p2 = obj.getPos2(t);
%             hold on
            g = plot([p1(1) p2(1)],[p1(2) p2(2)],'r-*',...
                'linewidth',2,'markersize',5);
%             plot([0 p1(1)],[0 p1(2)],'b-*',...
%                 'linewidth',2,'markersize',5)
%             hold off
        end
        
        function animation(obj)
            axis tight manual
            filename = 'plots/animation.gif';
            p1 = obj.pos1();
            p2 = obj.pos2();
            % Anfangskonfiguration zeichnen
            h_fig = figure(); clf;
            hold on
            plot([p1(1,1) p2(1,1)],[p1(2,1) p2(2,1)],'r-*',...
                'linewidth',2,'markersize',5)
            plot([0 p1(1,1)],[0 p1(2,1)],'b-*',...
                'linewidth',2,'markersize',5)
            hold off

            % Zeichenbereich festlegen
            axis equal
            %title('Ebener Manipulator')
            xRange = [-(obj.L1+obj.L2) (obj.L1+obj.L2)];
            yRange = xRange;
            %xRange = [-1 1];
            %yRange = [0 2];
            xlim(xRange)
            ylim(yRange)
            itm_formatfig(4,'FigSize',[10;10]);
            
            % Trajektorien plotten
%             hold on
%             plot(obj.z(1,:),obj.z(2,:),':b')
%             plot(obj.z(3,:),obj.z(4,:),':r')
            
            pause(0.5)

            % Animation-Schleife durchlaufen
            for i = 1:length(p1)
                % Capture the plot as an image 
                frame = getframe(h_fig); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                % Write to the GIF File 
                if i == 1 
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',obj.T0); 
                else 
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',obj.T0); 
                end 
                % Bild neu zeichnen
                figure(h_fig);cla; hold on
                hold on
                plot([p1(1,i) p2(1,i)],[p1(2,i) p2(2,i)],'r-*',...
                    'linewidth',2,'markersize',5)
                plot([0 p1(1,i)],[0 p1(2,i)],'b-*',...
                    'linewidth',2,'markersize',5)
                hold off
                %title(sprintf('Simulationszeit: %0.2fs', t_interp(idx)))

                % Zeichenbereich festhalten
                pause(obj.T0)
            end % for
        end
        
        function plot(obj, fig)
            %PLOT Plot the system states and control inputs
            %   fig: figure number
            if nargin < 2
                figure();
            else
                figure(fig);
            end
            %set(gcf, 'Position', get(0, 'Screensize'));
            % plot system states
            subplot(2,1,1);
            hold on
            title('Systemzustände');
            t = 0:(size(obj.x,2)-2);
            t = obj.T0 * t;
            xlabel('Zeit in s');
            yyaxis left
            plot(t,obj.x(1,1:end-1),'b-');
            plot(t,obj.x(2,1:end-1),'r-');
            ylabel('Winkel in rad');
            yyaxis right
            plot(t,obj.x(3,1:end-1),'b-.');
            plot(t,obj.x(4,1:end-1),'r-.');
            ylabel('Winkelgeschw. in rad/s');
            legend({'x1','x2','x3','x4'});%,'Location','south','NumColumns',2);
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'k';
            %ylim([-2 2]);
            %xlim([0 6000]);
            hold off
            % plot control inputs
            subplot(2,1,2);
            hold on
            title('Abweichung der Gesamtenergie');
            E = obj.energy(obj.x(:,1:end-1));
            plot(t,E-E(1))
            xlabel('Zeit in s');
            ylabel('Energiedifferenz in J');
            %ylim([-5 5]);
            %xlim([0 6000]);
            hold off
        end
        
        function plotStates(obj, fig)
            %PLOT Plot the system states and control inputs
            %   fig: figure number
            if nargin < 2
                figure();
            else
                figure(fig);
            end
            %set(gcf, 'Position', get(0, 'Screensize'));
            % plot system states
            hold on
            title('Systemzustände');
            t = 0:(size(obj.x,2)-2);
            t = obj.T0 * t;
            xlabel('Zeit in s');
            yyaxis left
            plot(t,obj.x(1,1:end-1),'b-');
            plot(t,obj.x(2,1:end-1),'r-');
            ylabel('Winkel in rad');
            yyaxis right
            plot(t,obj.x(3,1:end-1),'b-.');
            plot(t,obj.x(4,1:end-1),'r-.');
            ylabel('Winkelgeschw. in rad/s');
            legend({'x1','x2','x3','x4'});%,'Location','south','NumColumns',2);
            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'k';
            %ylim([-2 2]);
            %xlim([0 6000]);
            hold off
        end
        
        function plotEnergy(obj)
            t = 0:(size(obj.x,2)-2);
            t = obj.T0 * t;
            hold on
            %title('Abweichung der Gesamtenergie');
            E = obj.energy(obj.x(:,1:end-1));
            plot(t,E-E(1))
            xlabel('Zeit in s');
            ylabel('Energiedifferenz in J');
            %ylim([-5 5]);
            %xlim([0 6000]);
            hold off
        end
        
        function plotZ(obj, fig)
            %PLOT Plot the system states and control inputs
            %   fig: figure number
            if nargin < 2
                figure;
            else
                figure(fig);
            end
            %set(gcf, 'Position', get(0, 'Screensize'));
            % plot measurements
            c = linspace(1,10,size(obj.z,2));
            scatter(obj.z(1,:),obj.z(2,:),[],c);
        end
    end
end

