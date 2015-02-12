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
figure(5);
plot(XAxis,AmpX);

% Define tau steps
Size = 100;
tauStep = 1E-13;
tauAxis = 0:tauStep:(Size-1)*tauStep;
FStep = 1/tauAxis(end);
EgyStep = FStep*2*pi/meV2Hz;
EgyAxis = -EgyStep*(Size/2):EgyStep:EgyStep*(Size/2-1);
Egy_t = EgyAxis;
Egy_tau = flipud(EgyAxis');

% FRange = EgyRange*meV2Hz/(2*pi);% Total energy range
tauMax = tauAxis(end);     % Maximum scan time for tau

T = 2E-13;          % Delay T

% Full time axis parameters
TimeStart = -7E-13;

% Define matrix for final Y value
Y_All = zeros(Size,Size,3,XPts);
y0_1 = [0;0;0];

% Define tolerance for ode
options = odeset('RelTol',1e-4,'AbsTol',1e-8);
tic;

for j = 1:Size
    tau = tauAxis(j);
    Pos2 = Pos1 + tau;
    Pos3 = Pos2 + T;
    TimeMax = Pos3 + tauMax + 1E-12;
    PosMat = [Pos1;Pos2;Pos3];

    parfor k = 1:XPts
        Phase = 2*pi*XAxis(k).*k_pulse;
        Phase(1) = Phase(1)+pi;        % Give pi phase shift to puls A 
        XAmp = AmpX(k);
        y0 = cat(1,y0_1,Phase,PosMat,XAmp);
        [Time1,Y1] = ode23(@timeEvol_2Lvl,[TimeStart,TimeMax],y0,options);
        Y1 = Y1(:,1:3);
        
        % Take points after the 3rd pulse is incident
        Time1 = Time1 - Pos3;
        Pos3_Idx = find(Time1 > 0);
        Y1 = Y1(Pos3_Idx(1)-1:end,:);
        Time1 = Time1(Pos3_Idx(1)-1:end);
        
        Y = interp1(Time1',Y1,tauAxis');
        
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

dlmwrite('ESig_X1.dat',ESig_X,'\t');
size(ESig_X)
toc;
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
