function dy = timeEvol_2Lvl(t,y)

% y(4) = Pulse timing

% rho_11 = 1- rho_00
% y(1) = rho_11;
% In rotating frame, use sigma_01 instead of rho_01
% sigma_01 = y(2)+i*y(3)
% sigma_10 = y(2)-i*y(3)

global mu dt E BohrRad

meV2Hz = 241.79895E9*2*pi;      % hbar
e0 = 8.8541878E-12;     % epsilon0 (F/m)
hbar = 1.054572E-34;    % J s
% mu = 3.7E5;

% Define pulse parameters
% dt = 1.5E-13;       % Pulse width in s
% E = 5E-2*pi/mu;        % Pulse area

% System paramters
g01 = 0.73*meV2Hz;
G1 = 0.05*meV2Hz;
dw = 0*meV2Hz;

% Define MBE variables
EIS = 0*meV2Hz;
EID = 0*meV2Hz;
L = -(1/3)*mu*hbar/(e0*pi*BohrRad^3);         % Local field factor

% Define real and imaginary parts of E-field. Include amplitude along
% X-axis
E_R = E*(1/(dt*sqrt(pi)))*exp(-((t-y(4))/dt)^2);

E_I = 0;

H_Int = [-G1, -mu*E_I, -mu*E_R; mu*E_I, -g01, dw; mu*E_R, -dw, -g01];

H_MBE = [0, mu*L*y(3), -mu*L*y(2);...
    -mu*L*y(3), -EID*y(1), EIS*y(1);...
    mu*L*y(2), -EIS*y(1), -EID*y(1)];

HConst = [0; (1/2)*mu*(E_I - L*y(3)); (1/2)*mu*(E_R + L*y(2))];

dy = zeros(4,1);

dy(1:3) = (H_Int+H_MBE)*y(1:3) - HConst;