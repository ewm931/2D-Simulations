function dy = timeEvol_Inhom_2Lvl(t,y)

% y(end-6:end-4) = Phase
% y(end-3:end-1) = Pulse timing
% y(end) = E-field amplitude along X Axis

% rho_11 = 1- rho_00
% y(1) = rho_11;
% In rotating frame, use sigma_01 instead of rho_01
% sigma_01 = y(2)+i*y(3)
% sigma_10 = y(2)-i*y(3)

meV2Hz = 241.79895E9*2*pi;      % hbar

% Define pulse parameters
dt = 1.5E-13;       % Pulse width in s
E = 1E-1*pi;        % Pulse area

% System paramters
g01 = 0.1*meV2Hz;   % Dephasing rate
G1 = 0.05*meV2Hz;    % Population decay rate
dw_0 = 0*meV2Hz;    % Detuning of resonance center from laser center

% Many body terms
EIS = 0*meV2Hz;
EID = 0*meV2Hz;
L = 0.1*(1/(dt*sqrt(pi)));

% Define inhomogeneous width (Make sure same values are used in the main
% program)
g_Inhom = 0.4*meV2Hz;   % Inhomogeneity (sigma)
FPts = 25;
% CFId = 3*(FPts-1)/2+1;      % Center frequency rho_11 index
InhomRange = 8*g_Inhom;
dInhom = InhomRange/(FPts-1);
InhomEgy = -0.5*InhomRange : dInhom : 0.5*InhomRange;
% Set amplitude for each frequency component (as dipole factor)
Amp = exp(-0.5*(InhomEgy./g_Inhom).^2);
Amp = Amp./sum(Amp);

dw = InhomEgy + dw_0;

% Define real and imaginary parts of E-field. Include amplitude along
% X-axis
Idx = 3*FPts+1;     % Define start index for pulse parameters
E_R = (1/2)*y(Idx+6)*E*(1/(dt*sqrt(pi)))*(exp(-((t-y(Idx+3))/dt)^2).*cos(y(Idx)) +...
    exp(-((t-y(Idx+4))/dt)^2).*cos(y(Idx+1)) + exp(-((t-y(Idx+5))/dt)^2).*cos(y(Idx+2)));

E_I = (1/2)*y(Idx+6)*E*(1/(dt*sqrt(pi)))*(exp(-((t-y(Idx+3))/dt)^2).*sin(y(Idx)) +...
    exp(-((t-y(Idx+4))/dt)^2).*sin(y(Idx+1)) + exp(-((t-y(Idx+5))/dt)^2).*sin(y(Idx+2)));

% Interaction hamiltonian w/ only E-field
HConst1 = [0; E_I; E_R];
HConst = repmat(HConst1,FPts,1);
% Interaction hamiltonian matrix including MBE terms
H_Int = [-G1, -2*E_I, -2*E_R; 2*E_I, -g01, dw(1); 2*E_R, -dw(1), -g01];
H_MBE = [0, 2*Amp(1)*L*y(3), -2*Amp(1)*L*y(2);...
    -2*Amp(1)*L*y(3), -Amp(1)*EID*y(1), Amp(1)*EIS*y(1);...
    2*Amp(1)*L*y(2), -Amp(1)*EIS*y(1), -Amp(1)*EID*y(1)];
% LF term with only E-field part of the Hamiltonian
HConst_LF = [0; -Amp(1)*L*y(3); Amp(1)*L*y(2)];
for j = 2:FPts
    H_Int1 = [-G1, -2*E_I, -2*E_R; 2*E_I, -g01, dw(j); 2*E_R, -dw(j), -g01];
    H_MBE1 = [0, 2*Amp(j)*L*y(3*(j-1)+3), -2*Amp(j)*L*y(3*(j-1)+2);...
        -2*Amp(j)*L*y(3*(j-1)+3), -Amp(j)*EID*y(3*(j-1)+1), Amp(j)*EIS*y(3*(j-1)+1);...
        2*Amp(j)*L*y(3*(j-1)+2), -Amp(j)*EIS*y(3*(j-1)+1), -Amp(j)*EID*y(3*(j-1)+1)];
    HConst_LF1 = [0; -Amp(j)*L*y(3*(j-1)+3); Amp(j)*L*y(3*(j-1)+2)];
    H_Int = blkdiag(H_Int,H_Int1);
    H_MBE = blkdiag(H_MBE,H_MBE1);
    HConst_LF = cat(1,HConst_LF,HConst_LF1);
end
HConst_Tot = HConst + HConst_LF;

dy = zeros(FPts*3+7,1);

dy(1:3*FPts) = (H_Int + H_MBE)*y(1:3*FPts) - HConst_Tot;