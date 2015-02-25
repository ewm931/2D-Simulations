%% Numerically simulate 2D spcectrum
% Simple 2-level system. Phase matching is implemented.
% Changing the pulse duration, time step and total time can be tricky.
% Currently I'm using high dephasing rate and population time to keep the
% simulation time small.

%% Method
% Define 3 Gaussian pulses. Repeat for different tau values.
% Since tau changes, the total time duration changes.
% Simulating only S1 spectrum. Hence, will keep T fixed.
% ODEs for different x positions are solved in a loop. The phase term for a
% specific x-position is passed as a variable to enable parallel execution
% of the loop.

clear all;% clc;

global mu dt E BohrRad

meV2Hz = 241.79895E9*2*pi;      % 1/hbar
mu = 1.4E6;           % Actually mu/hbar
dt = 1.5E-13;       % Pulse width in s
E = 5E-2*pi/mu;        % Peak electric field amplitude

% Normalization constants (all in SI units)
hbar = 1.054572E-34;    % J s
BohrRad = 8.3E-9;       % 10 nm
Res_w = 2*pi*3.75E14;   % Resonance angular frequency
c = 2.9979246E8;        % m/s
e0 = 8.8541878E-12;     % epsilon0 (F/m)

%% Solve OBE's during tau

Pos = 0;           % Timing of the pulse

% Define tau steps
Size = 100;
tStep = 1E-13;
tAxis = 0:tStep:(Size-1)*tStep;
FStep = 1/tAxis(end);
EgyStep = FStep*2*pi/meV2Hz;
EgyAxis = -EgyStep*(Size/2):EgyStep:EgyStep*(Size/2-1);

tPulse = tAxis(Size/2+1);
E_R = E*(1/(dt*sqrt(pi)))*exp(-((tAxis-tPulse)/dt).^2);
E_R = ifftshift(E_R);
E_R_w = fftshift(fft(E_R));
% figure(100);
% plot(tAxis,E_R);
% FRange = EgyRange*meV2Hz/(2*pi);% Total energy range
tMax = tAxis(end);     % Maximum scan time for tau

% Full time axis parameters
TimeStart = -7E-13;

% Define matrix for final Y value
y0_1 = [0;0;0];

% Define tolerance for ode
options = odeset('RelTol',1e-6,'AbsTol',1e-12);

TimeMax = Pos + tMax + 1E-12;
PosMat = Pos;

y0 = cat(1,y0_1,PosMat);
[Time1,Y1] = ode23(@timeEvol_2Lvl,[TimeStart,TimeMax],y0,options);
Y1 = Y1(:,1:3);

% Take points after the 3rd pulse is incident
Time1 = Time1 - Pos;
Pos_Idx = find(Time1 > 0);
Y1 = Y1(Pos_Idx(1)-1:end,:);
Time1 = Time1(Pos_Idx(1)-1:end);

Y = interp1(Time1',Y1,tAxis');

figure(1);
subplot(131);
plot(tAxis,Y(:,1));
title('Population');
subplot(132);
plot(tAxis,Y(:,2));
title('Re(\rho_{01})');
subplot(133);
plot(tAxis,Y(:,3));
title('Im(\rho_{01})');

%%
% clear Y_All
%%
% Save 3rd order polarization (sigma_10)
Pol = Y(:,2) - 1i*Y(:,3);

ESig = 1i*hbar*mu*Res_w*(1/(2*c*e0*pi*BohrRad^2))*Pol;

ESig_w = fftshift(fft(ESig));

% E_Tot = E_R' + ESig;
E_Tot_w =  E_R_w' + ESig_w;
Trans_Int = abs(E_Tot_w).^2;
Inc_Int = abs(E_R_w').^2;
Abs_Int = Inc_Int - Trans_Int;
figure(3);
plot(EgyAxis,Abs_Int);
figure(4);
plot(EgyAxis,Inc_Int,EgyAxis,Trans_Int);

% figure(5);
% plot(tAxis,ESig,tAxis,

