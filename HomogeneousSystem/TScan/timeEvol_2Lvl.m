function dy = timeEvol_2Lvl(t,y)

% y(4:6) = Phase
% y(7:9) = Pulse timing
% y(10) = E-field amplitude along X Axis

% rho_11 = 1- rho_00
% y(1) = rho_11;
% In rotating frame, use sigma_01 instead of rho_01
% sigma_01 = y(2)+i*y(3)
% sigma_10 = y(2)-i*y(3)

meV2Hz = 241.79895E9*2*pi;      % hbar

% Define pulse parameters
dt = 1.5E-13;       % Pulse width in s
E = 5E-2*pi;        % Pulse area

% System paramters
g01 = 0.4*meV2Hz;
G1 = 0.2*meV2Hz;
dw = 0*meV2Hz;

% Define MBE variables
EIS = 0*meV2Hz;
EID = 0*meV2Hz;
L = -0.5*(1/(dt*sqrt(pi)));         % Local field factor

% Define real and imaginary parts of E-field. Include amplitude along
% X-axis
E_R = (1/2)*y(10)*E*(1/(dt*sqrt(pi)))*(exp(-((t-y(7))/dt)^2).*cos(y(4)) +...
    exp(-((t-y(8))/dt)^2).*cos(y(5)) + exp(-((t-y(9))/dt)^2).*cos(y(6)));

E_I = (1/2)*y(10)*E*(1/(dt*sqrt(pi)))*(exp(-((t-y(7))/dt)^2).*sin(y(4)) +...
    exp(-((t-y(8))/dt)^2).*sin(y(5)) + exp(-((t-y(9))/dt)^2).*sin(y(6)));

H_Int = [-G1, -2*E_I, -2*E_R; 2*E_I, -g01, dw; 2*E_R, -dw, -g01];

% H_MBE = [0, 0, 0; 0, -EID*y(1), EIS*y(1);...
%     0, -EIS*y(1), -EID*y(1)];

H_MBE = [0, 2*L*y(3), -2*L*y(2);...
    -2*L*y(3), -EID*y(1), EIS*y(1);...
    2*L*y(2), -EIS*y(1), -EID*y(1)];

HConst = [0; E_I - L*y(3); E_R + L*y(2)];

dy = zeros(10,1);

dy(1:3) = (H_Int+H_MBE)*y(1:3) - HConst;