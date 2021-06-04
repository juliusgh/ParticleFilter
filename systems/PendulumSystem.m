classdef PendulumSystem < System
    %PENDULUMSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    properties
        T % sampling time
        G % gravity constant
        L % pendulum length
        m % mass
        J % Trägheitsmoment
    end
    
    methods
        function obj = PendulumSystem()
            %SYSTEM Construct an instance of this class
            %   Simple mechanical test system
            obj.dim = 4; % system dimensions
            obj.dimz = 2; % measurement dimensions
            obj.dimu = 2; % control input dimensions
            obj.x0 = [pi/2;0;0;0]; % initial condition
            obj.V = HyperRect([-0.01 0.01],[-0.01 0.01],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.W = HyperRect([-0.3 0.3],[-0.3 0.3],[-0.3 0.3],[-0.3 0.3]); % measurement noise support
            obj.Vsim = HyperRect([-0.01 0.01],[-0.01 0.01],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.Wsim = HyperRect([-0.3 0.3],[-0.3 0.3],[-0.3 0.3],[-0.3 0.3]); % measurement noise support
            obj.X0 = obj.x0 + obj.W; % initial noise support
            obj.uk = @(k) [0;0]; % fixed input function
            obj.T = 0.01;
            obj.G = 9.81;
            obj.L = 1;
            obj.m = 1;
            obj.J = obj.m * obj.L^2;
        end
        
        function x_ = g(obj,x,u)
            %G process function with input
            %   x: current system state, u: control input
%             x_ = obj.gpredict(x,u);
            if size(x,2) > 1
                x_ = NaN(size(x));
                for i = 1:size(x,2)
                    x_(:,i) = obj.g(x(:,i),u);
                end
            else
                options = odeset('MaxStep',1e-3);
                [~,y] = ode45(@(t,x) obj.f(t,x,u),[0 obj.T],x,options);
                x_ = y(end,:)';
            end
        end
        
        function x_ = gpredict(obj,x,u)
            %G process function with input for prediction (less accurate)
            %   x: current system state, u: control input
            maxStep = 0.001;
            n = ceil(obj.T/maxStep) - 1;
            for i=1:n
                x = obj.rk45(@(t,x) obj.f(t,x,u),0,maxStep,x);
            end
            x_ = obj.rk45(@(t,x) obj.f(t,x,u),0,obj.T-n*maxStep,x);
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
            [x1,x2,x3,x4] = unpack(x);
            z = zeros(size(x1));
            upart = [zeros(2);eye(2)]*u;
            if isa(x,'interval')
                upart = upart';
            end
            x_ = [x3;x4;z;z] + obj.n(x) + upart;
        end
        
        function n = n(obj,x)
            %N non-linear process function part
            %   x: current system state
            [x1,x2,x3,x4] = unpack(x);
            z = zeros(size(x1));
            s1 = sin(x1);
            s2 = sin(x2);
            s12 = sin(x1-x2);
            c12 = cos(x1-x2);
            n = 1./(2-c12.^2) .* (s12.*[z;z;-x4.^2-x3.^2.*c12;2*x3.^2+x4.^2.*c12]...
                +obj.G/obj.L.*[z;z;s2.*c12-2.*s1;2.*s1.*c12-2.*s2]);
        end
        
        function z = h(~,x)
            %H measurement function
            %   x: current system state
            %z = [eye(2) zeros(2)]*x;
            z = x;
        end
        
        function u = uref(obj,xref)
            %UREF get desired control input for reference state
            %   x: reference system state
            xr = obj.f(0,xref,[0;0]);
            u = -xr(3:4,:);
        end
        
        function E = energy(obj,x)
            %ENERGY calculate total energy of the system
            %   x: system states
            if nargin < 2
                x = obj.x;
            end
            Ekin = 0.5 * obj.m * obj.L^2 * (2*x(3,:).^2 + x(4,:).^2 + 2*x(3,:).*x(4,:).*cos(x(1,:)-x(2,:)));
            Epot = obj.m * obj.G * obj.L * (4 - 2*cos(x(1,:)) - cos(x(2,:)));
            E = Ekin + Epot;
        end
        
        function p = pos1(obj,x)
            if nargin < 2
                x = obj.x;
            end
            p = [obj.L * sin(x(1,:)); -obj.L * cos(x(1,:))];
        end
        
        function p = pos2(obj,x)
            if nargin < 2
                x = obj.x;
            end
            p = [obj.L*sin(x(1,:))+obj.L*sin(x(2,:));-obj.L*cos(x(1,:))-obj.L*cos(x(2,:))];
        end
        
        function animation(obj)
            p1 = obj.pos1();
            p2 = obj.pos2();
            % Anfangskonfiguration zeichnen
            h_fig = figure(2); clf;
            hold on
            plot([p1(1,1) p2(1,1)],[p1(2,1) p2(2,1)],'r-*',...
                'linewidth',2,'markersize',5)
            plot([0 p1(1,1)],[0 p1(2,1)],'b-*',...
                'linewidth',2,'markersize',5)
            hold off

            % Zeichenbereich festlegen
            axis equal
            title('Ebener Manipulator')
            xRange = [-2*obj.L 2*obj.L];
            yRange = xRange;
            xlim(xRange)
            ylim(yRange)
            
            % Trajektorien plotten
            hold on
            plot(obj.z(1,:),obj.z(2,:),':b')
            plot(obj.z(3,:),obj.z(4,:),':r')
            
            pause(0.5)

            % Animation-Schleife durchlaufen
            for i = 1:length(p1)
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
                pause(obj.T)
            end % for
        end
        
        function x_ = f2(obj,x)
            s1 = sin(x(1,:));
            s2 = sin(x(2,:));
            s12 = sin(x(1,:)-x(2,:));
            c12 = cos(x(1,:)-x(2,:));
            x_ = [2 c12;c12 1]\(s12*[-x(4)^2;x(3)^2]-obj.G/obj.L*[2*s1;s2]);
        end
        
        function x_ = f1(obj,a1,a2,a11,a22,a111,a222)
            x_ = obj.L*[2 cos(a1-a2);cos(a1-a2) 1]*[a111;a222]+obj.L*[a22^2*sin(a1-a2);-a11^2*sin(a1-a2)]+obj.G*[2*sin(a1);sin(a2)];
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
            t = obj.T * t;
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
            plot(t,obj.energy(obj.x(:,1:end-1)))
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

