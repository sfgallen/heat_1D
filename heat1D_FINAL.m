% model numerically solves the transient conduction-advection heat equation
% Equations are after Moore and England (2001) and Braun, van der Beek, and
% Batt (2006).
% Author: Sean F. Gallen
% email: sean.gallen[at]colostate.edu
% Date modified: 12/07/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set parameters
run_time = 250;     % run time in Ma.
ei = 0.5;           % initial erosion rate in km/Myr
Ts = 10;            % surface temperature in degrees C
Tb = 1300;          % temperature at the base of the lithosphere
k = 25;             % thermal conduction coefficient
lambda = 0.4;       % this is for stability, usually between 0.5 and 0.05 is good.
                    % reduce the number of the problem becomes unstable
                    
Zmin = 0;           % min Z in terms of km
Zmax = 100;         % max Z in terms of km (approximation of base of lithosphere
dZ = 0.25;          % increments of Z in km -- as this decreases lambda needs to decrease

n_plots = 25;       % number of transient thermal geotherms plotted

%%
% derive parameters for lithospheric heat production profile after Hasterok 
% and Chapman (2006)
K = 25;             % thermal conductivity
P = 0.6;            % Partition Coefficient Pollack and Chapman (1977)
qo = 75;            % surface hear flow
qr = P*qo;          % reduced heat flow
D = 8;              % charateristic thickness of heat producing layer
Zc = 35;            % crustal thickness
Zlc = 17.5;         % depth to lower crust
Alc = 0.4;          % heat production in lower crust
Am = 0.02;          % heat production in mantle
Ao = (1-P)*qo/D;    % surface heat production

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = Zmin:dZ:Zmax;   % Z (depth) vector

% derive heat production profile 
A = Ao*exp(-Z/D);
A(Z >= Zlc) = Alc;
A(Z >= Zc) = Am;
AK = A./K;          % set this to zero if you want to ignore lithospheric heat production

%%
% estimate To with ei or ~0 with Brauns peclet number equation
Pe = (1e-9.*max(Z))/k; % peclect number from braun
To = (1-exp(1).^(-Pe.*(Z./max(Z))))./(1-exp(1).^(-Pe));
To = To.*(Tb-Ts)+Ts;

% make dt such that lambda is 0.5
dt = lambda*(dZ^2)/k;

% cast time vector for numerical calculations
time_v = 0:dt:run_time;

figure(1)
hold off
plot(To,Z,'k-'); hold on
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
    dTdt = ddTddz+dTdz+(AK);
    
    % add to previous thermal state
    To = To + dTdt;

    % update the finite difference vector
    To_pad(2:end-1) = To;
    
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
