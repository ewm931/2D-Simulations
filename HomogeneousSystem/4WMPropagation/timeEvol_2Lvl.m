function dy = timeEvol_2Lvl(t,y,tAxis,ESig_X,mu,Lfact,EISfact,EIDfact,gamma,E,probe)
% mu, L, EIS, EID, gamma, E = pulse area

% y(4:6) = Phase
% y(7:9) = Pulse timing
% y(10) = E-field amplitude along X Axis
% y(11) = Real E-field due to previous wells
% y(12) = Imag E-field due to previous wells

% rho_11 = 1- rho_00
% y(1) = rho_11;
% In rotating frame, use sigma_01 instead of rho_01
% sigma_01 = y(2)+i*y(3)
% sigma_10 = y(2)-i*y(3)

meV2Hz = 241.79895E9*2*pi;      % hbar

% Define pulse parameters
dt = 1.5E-13;       % Pulse width in s
%E = 5E-2*pi;        % Pulse area

% System paramters
g01 = gamma*meV2Hz; %0.5
G1 = 0.05*meV2Hz;
dw = 0*meV2Hz;

% Define MBE variables
EIS = EISfact*meV2Hz;
EID = EIDfact*meV2Hz;
L = mu*Lfact;         % Local field factor

% Intepolate ESig_X(t_Axis) to get Ep_R(t) and Ep_I(t)
Ep_R = interp1(tAxis,real(ESig_X),t,'linear',0);
Ep_I = interp1(tAxis,imag(ESig_X),t,'linear',0);

% Define real and imaginary parts of E-field. Include amplitude along
% X-axis
E_R = y(10)*E*(1/(dt*sqrt(pi)))*(1/mu)*(exp(-((t-y(7))/dt)^2).*cos(y(4)) +...
    exp(-((t-y(8))/dt)^2).*cos(y(5)) + probe*exp(-((t-y(9))/dt)^2).*cos(y(6))) + Ep_R;

E_I = y(10)*E*(1/(dt*sqrt(pi)))*(1/mu)*(exp(-((t-y(7))/dt)^2).*sin(y(4)) +...
    exp(-((t-y(8))/dt)^2).*sin(y(5)) + probe*exp(-((t-y(9))/dt)^2).*sin(y(6))) + Ep_I;

H_Int = [-G1, -mu*E_I, -mu*E_R; mu*E_I, -g01, dw; mu*E_R, -dw, -g01];

% H_MBE = [0, 0, 0; 0, -EID*y(1), EIS*y(1);...
%     0, -EIS*y(1), -EID*y(1)];

H_MBE = [0, mu*L*y(3), -mu*L*y(2);...
    -mu*L*y(3), -EID*y(1), EIS*y(1);...
    mu*L*y(2), -EIS*y(1), -EID*y(1)];

HConst = [0; (1/2)*mu*(E_I - L*y(3)); (1/2)*mu*(E_R + L*y(2))];

dy = zeros(10,1);

dy(1:3) = (H_Int+H_MBE)*y(1:3) - HConst;