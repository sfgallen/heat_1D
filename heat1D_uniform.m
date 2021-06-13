% model numerically solves the transient conduction-advection heat equation
% Equations are after Moore and England (2001) and Braun, van der Beek, and
% Batt (2006).
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% Date modified: 12/07/2019

%% clear work space
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set parameters
run_time = 100;     % run time in Ma.
ei = 0.5;           % initial erosion rate in km/Myr
Ts = 10;            % surface temperature in degrees C
Tb = 540;           % temperature at the base of the crust
k = 25;             % thermal conduction coefficient
A = 4.5;            % crustal heat production in C/Myr
lambda = 0.4;       % this is for stability, usually between 0.5 and 0.05 is good.
                    % reduce the number of the problem becomes unstable
                   
Zmin = 0;           % min Z in terms of km
Zmax = 50;         % max Z in terms of km (approximation of base of lithosphere
dZ = 0.25;          % increments of Z in km -- as this decreases lambda needs to decrease

Z = Zmin:dZ:Zmax;   % Z (depth) vector

n_plots = 25;       % number of transient thermal geotherms plotted

%%
% estimate To with ei or ~0 with Brauns peclet number equation
Pe = (1e-9.*max(Z))/k; % peclect number from braun
To = (1-exp(1).^(-Pe.*(Z./max(Z))))./(1-exp(1).^(-Pe));
To = To.*(Tb-Ts)+Ts;

% Batt Brandon (2002) final steady-state solution
Tz = Ts + (Tb-Ts+((A.*Zmax)./ei)).*((1-exp(-ei.*Z./k))./(1-exp(-ei.*Zmax./k)))-((A.*Z)./ei);

% make dt such that lambda is 0.5
dt = lambda*(dZ^2)/k;

% cast time vector for numerical calculations
time_v = 0:dt:run_time;

figure(1)
hold off
plot(To,Z,'k-'); hold on
plot(Tz,Z,'k--','linewidth',2,'color',[0.5 0.5 0.5]);
axis ij
xlabel('Temperature (^oC)')
ylabel('Depth (km)')

% set up plotting for transient profiles
t_plots = round(round(length(time_v))./n_plots);
n_plots = round(round(length(time_v))./t_plots);
np = 1;

% colors for the transient profiles
cols = hsv(n_plots);

% place thermal boundary conditions on To for curvature calculation and
% slope calculation. This will take care of edge effects.
To_pad = nan(1,length(To)+2);
To_pad(2:end-1) = To;
To_pad(1) = Ts;
To_pad(end) = Tb;

h = waitbar(0,'running thermal code...');
for t = 1:length(time_v)

    % calculate the first derivative and advect;
    dTdz = ((To_pad(2:end-1) - To_pad(3:end))./dZ).*-ei.*dt;
    
    % calculate second derivative and diffuse
    ddTddz = ((To_pad(3:end) - 2.*(To_pad(2:end-1)) + To_pad(1:end-2))).*lambda;
    
    % add the conduction and advection to gether to predict change in
    % thermal profile
    dTdt = ddTddz+dTdz+(A.*dt);
    
    % add to previous thermal state
    To = To + dTdt;

    % update the finite difference vector and deal with boundaries
    To_pad(2:end-1) = To;
    To(1) = Ts;
    To(end) = Tb;
    
    % plot as needed
    if rem(t,t_plots) == 0
        plot(To,Z,'-','color',cols(np,:));
        np = np+1;
    end
    waitbar(t/length(time_v),h);
end
close(h)

% final model position
plot(To,Z,'k-','linewidth',2); hold on

