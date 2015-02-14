function [] = TwoLevel_2D_Fn(Index)
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

% clear all;% clc;

meV2Hz = 241.79895E9*2*pi;      % hbar

%% Solve OBE's during tau

Pos1 = 0;           % Timing of 1st pulse

% Define X-axis for taking FT to implement phase matching
XPts = 101;          % # of points along X axis
dx = 0.05;
XAxis1 = -0.5*dx*(XPts-1):dx:0.5*dx*(XPts-1);
XAxis = ifftshift(XAxis1);
k_pulse = [-4;1;-1];

% Decay in X
wX = 1;
AmpX = exp(-(XAxis./wX).^2);
% figure(5);
% plot(XAxis,AmpX);

% Define T and t axes
Size_t = 20;
tStep = 1E-13;
tAxis = 0:tStep:(Size_t-1)*tStep;
Size_T = 10;
TStep = 1E-13;
TAxis = 0:TStep:(Size_T-1)*TStep;

switch Index
    case 0
        PhaseA = 0;
        PhaseB = 0;
        
        % Define k axis
        dk = 1/(max(XAxis)-min(XAxis));
        kMax = 1/dx;
        kAxis = 0:dk:kMax;
        dlmwrite('kAxis.dat',kAxis,'\t');
        dlmwrite('WVect.dat',k_pulse,'\t');
        
        % Define frequency/energy axes
        F_tStep = 1/tAxis(end);
        Egy_tStep = F_tStep*2*pi/meV2Hz;
        Egy_tAxis = -Egy_tStep*(Size_t/2):Egy_tStep:Egy_tStep*(Size_t/2-1);
        Egy_t = Egy_tAxis;
        
        F_TStep = 1/TAxis(end);
        Egy_TStep = F_TStep*2*pi/meV2Hz;
        Egy_TAxis = -Egy_TStep*(Size_T/2):Egy_TStep:Egy_TStep*(Size_T/2-1);
        Egy_T = Egy_TAxis;
        
        dlmwrite('Egy_t.dat',Egy_t,'\t');
        dlmwrite('Egy_capT.dat',Egy_T,'\t');
        dlmwrite('tAxis.dat',tAxis,'\t');
        dlmwrite('capTAxis.dat',TAxis,'\t');
        
    case 1
        PhaseA = pi;
        PhaseB = 0;
        
    case 2
        PhaseA = pi;
        PhaseB = pi;
        
    case 3
        PhaseA = 0;
        PhaseB = pi;
end

tMax = tAxis(end);     % Maximum scan time for tau

tau = 2E-13;          % Delay tau

% Full time axis parameters
TimeStart = -7E-13;

% Define matrix for final Y value
NVar = 3;
Y_All = zeros(Size_T,Size_t,NVar,XPts);
y0_1 = zeros(NVar,1);

% Define tolerance for ode
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
tic;

for j = 1:Size_T
    T = TAxis(j);
    Pos2 = Pos1 + tau;
    Pos3 = Pos2 + T;
    TimeMax = Pos3 + tMax + 1E-12;
    PosMat = [Pos1;Pos2;Pos3];

    parfor k = 1:XPts
        Phase = 2*pi*XAxis(k).*k_pulse;
        Phase(1) = Phase(1) + PhaseA;
        Phase(2) = Phase(2) + PhaseB;
        XAmp = AmpX(k);
        y0 = cat(1,y0_1,Phase,PosMat,XAmp);
        [Time1,Y1] = ode23(@timeEvol_2Lvl,[TimeStart,TimeMax],y0,options);
        Y1 = Y1(:,1:NVar);
        
        % Take points after the 3rd pulse is incident
        Time1 = Time1 - Pos3;
        Pos3_Idx = find(Time1 > 0);
        Y1 = Y1(Pos3_Idx(1)-1:end,:);
        Time1 = Time1(Pos3_Idx(1)-1:end);
        
        Y = interp1(Time1',Y1,tAxis');
        
        Y_All(j,:,:,k) = Y;
        
    end
    j
%     length(Time1)
end
toc;
%%
% clear Y_All
%%
% Save 3rd order polarization (sigma_10)
Pol3_X = Y_All(:,:,2,:) - 1i*Y_All(:,:,3,:);

ESig_X = 1i*Pol3_X;

OutFile = strcat('ESig_X',num2str(Index),'.dat');
dlmwrite(OutFile,ESig_X,'\t');
% size(ESig_X)
toc;
end
% ESig_K = fft(ESig_X,[],4);
% proj1 = sum(sum(abs(ESig_K),2),1);
% proj1 = reshape(proj1,1,[]);
% figure(1);
% semilogy(kAxis,proj1);
% 
% %%
% % Select out the 3rd order signal
% kSig = -k_pulse(1)+k_pulse(2)+k_pulse(3);
% % kSig = 3.1;
% kDiff = abs(kAxis - kSig);
% L = find(kDiff == min(kDiff));
% 
% % % Cut relevant points in k-space
% % NkPts = 2;
% ESig_K_Cut = zeros(size(ESig_K));
% % ESig_K_Cut(:,:,:,L-NkPts:L+NkPts) = ESig_K(:,:,:,L-NkPts:L+NkPts);
% 
% % Apply window to smooth the function
% WinSize = 5;            % # points on one side of peak
% % Win = hamming(2*WinSize+1);
% % Win = ones(1,2*WinSize+1);
% Window = zeros(1,XPts);
% Window(L-WinSize:L+WinSize) = 1;
% parfor j = 1:XPts
%     ESig_K_Cut(:,:,:,j) = ESig_K(:,:,:,j).*Window(j);
% end
% 
% proj1 = sum(sum(abs(ESig_K_Cut),2),1);
% proj1 = reshape(proj1,1,[]);
% figure(2);
% plot(kAxis,proj1);
% 
% 
% % IFT to spatial domain
% ESig_X_Cut = ifft(ESig_K_Cut,[],4);
% ESig_t = ESig_X_Cut(:,:,:,1);
% ESig_t = reshape(ESig_t,Size,Size);
% % ESig_t = interp2(tAxis,tauAxis',ESig_t,tAxis,tAxis');
% %% Define plot ranges
% TimeMax = 1E-11;    % seconds
% EgyMax = 7;         % meV (one side)
% TimeRange = [0 TimeMax 0 TimeMax];
% EgyRange = [-EgyMax EgyMax -EgyMax EgyMax];
% NCont = 50;
% 
% VMax = max(max(abs(ESig_t)));
% figure(3)
% set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 1 4 14 5 ]);
% subplot(131);
% contourf(tauAxis,tauAxis',abs(ESig_t),linspace(0,VMax,NCont),'LineStyle','none');
% xlabel('t (s)');
% ylabel('\tau (s)');
% title('Absolute value');
% axis(TimeRange);
% subplot(132);
% contourf(tauAxis,tauAxis',real(ESig_t),linspace(-VMax,VMax,NCont),'LineStyle','none');
% xlabel('t (s)');
% ylabel('\tau (s)');
% title('Real part');
% axis(TimeRange);
% subplot(133);
% contourf(tauAxis,tauAxis',imag(ESig_t),linspace(-VMax,VMax,NCont),'LineStyle','none');
% xlabel('t (s)');
% ylabel('\tau (s)');
% title('Imaginary part');
% axis(TimeRange);
% 
% % ESig_t = ESig_t';
% % NFin = 512;
% % tauMax1 = 1E-11;     % Maximum scan time for tau
% % tauStep1 = tauMax1/NFin;
% % tauAxis1 = 0:tauStep1:tauMax1-tauStep1;
% % tAxis1 = tauAxis1';
% % ESig_t1 = interp2(tauAxis1,tAxis,ESig_t,tauAxis1,tAxis1);
% % Evaluate frequency-domain signal
% Spec_2D = fftshift(fft(fftshift(fft(ESig_t,[],2),2),[],1),1);
% Spec_2D = flipud(Spec_2D);
% % Spec_2D = abs(Spec_2D);
% VMax1 = max(max(abs(Spec_2D)));
% % Ft = FAxis - FAxis(Size/2+1);
% % Egy_t = EgyAxis;
% % Ftau = -(FAxis - FAxis(Size/2))';
% % Egy_tau = flipud(EgyAxis');
% figure(4);
% set(gcf, 'Units', 'inch');
% set(gcf, 'position', [ 1 4 12 5 ]);
% subplot(131);
% contourf(Egy_t,Egy_tau,abs(Spec_2D),linspace(0,VMax1,NCont),'LineStyle','none');
% xlabel('\omega_t (meV)');
% ylabel('\omega_\tau (meV)');
% title('Absolute value');
% axis(EgyRange);
% subplot(132);
% contourf(Egy_t,Egy_tau,real(Spec_2D),linspace(-VMax1,VMax1,NCont),'LineStyle','none');
% xlabel('\omega_t (meV)');
% ylabel('\omega_\tau (meV)');
% title('Real part');
% axis(EgyRange);
% subplot(133);
% contourf(Egy_t,Egy_tau,imag(Spec_2D),linspace(-VMax1,VMax1,NCont),'LineStyle','none');
% xlabel('\omega_t (meV)');
% ylabel('\omega_\tau (meV)');
% title('Imaginary part');
% axis(EgyRange);
