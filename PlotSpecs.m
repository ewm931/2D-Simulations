% Define date and run index
% Date = '2014_12_25';
% Index = 'Run2';
% Folder = strcat('./',Date,'/',Index,'/');
global Folder
meV2Hz = 241.79895E9*2*pi;      % hbar

% Read axis files
kAxis = dlmread(strcat(Folder,'kAxis.dat'),'\t');
k_pulse = dlmread(strcat(Folder,'WVect.dat'),'\t');
tauAxis1 = dlmread(strcat(Folder,'tauAxis.dat'),'\t');

Size1 = size(tauAxis1,2);    % # of time points from simulation
XPts = size(kAxis,2);

%% Set time and energy axes after zero padding
Size = 1024;        % Set number points after zero padding
tauStep = tauAxis1(2)-tauAxis1(1);
tauAxis = 0:tauStep:(Size-1)*tauStep;
FStep = 1/tauAxis(end);
EgyStep = FStep*2*pi/meV2Hz;
EgyAxis = -EgyStep*(Size/2):EgyStep:EgyStep*(Size/2-1);
Egy_t = EgyAxis;
Egy_tau = flipud(EgyAxis');
dlmwrite(strcat(Folder,'Egy_t_Interp.dat'),Egy_t,'\t');
dlmwrite(strcat(Folder,'Egy_tau_Interp.dat'),Egy_tau,'\t');

tauAxis = tauAxis.*1E12;    % Convert to ps

%% Read files
ESig_X = dlmread(strcat(Folder,'ESig_X.dat'),'\t');
ESig_X = reshape(ESig_X,Size1,Size1,[]);
ESig_K = fft(ESig_X,[],3);

proj1 = sum(sum(abs(ESig_K),2),1);
proj1 = reshape(proj1,1,[]);
figure(1);
semilogy(kAxis,proj1);
xlabel('k');
ylabel('Polarization amplitude');
saveas(gcf, strcat(Folder, 'WaveVect'), 'emf');

% %% Define finer energy mesh
% Size = 1000;
% EgyRange = max(Egy_t1)-min(Egy_t1);
% EgyStep = EgyRange/Size;
% Egy_t = -Size/2*EgyStep:EgyStep:(Size/2-1)*EgyStep;
% Egy_tau = Egy_t';

%% K-time 2Dplots
ESig_K_t = sum(abs(ESig_K),1);
ESig_K_t = reshape(ESig_K_t,Size1,XPts);
ESig_K_tau = sum(abs(ESig_K),2);
ESig_K_tau = reshape(ESig_K_tau,Size1,XPts);
VMax1 = max(max(ESig_K_t));
VMax2 = max(max(ESig_K_tau));
NCont = 50;
figure(6);
set(gcf, 'Units', 'inch');
set(gcf, 'position', [ 0.5 1 10 5 ]);
subplot(121);
contourf(kAxis,tauAxis1',ESig_K_t,linspace(0,VMax1,NCont),'LineStyle','none');
xlabel('k');
ylabel('t (ps)');
title('k-t Projection');
colorbar ('Location','NorthOutside');
subplot(122);
contourf(kAxis,tauAxis1',ESig_K_tau,linspace(0,VMax2,NCont),'LineStyle','none');
xlabel('k');
ylabel('\tau (ps)');
title('k-\tau Projection');
colorbar ('Location','NorthOutside');
saveas(gcf, strcat(Folder, 'K-time'), 'emf');

%%
% Select out the 3rd order signal
kSig = -k_pulse(1)+k_pulse(2)+k_pulse(3);
% kSig = 3.1;
kDiff = abs(kAxis - kSig);
L = find(kDiff == min(kDiff));

% % Cut relevant points in k-space
% NkPts = 2;
ESig_K_Cut = zeros(size(ESig_K));
% ESig_K_Cut(:,:,:,L-NkPts:L+NkPts) = ESig_K(:,:,:,L-NkPts:L+NkPts);

% Apply window to smooth the function
WinSize = 5;            % # points on one side of peak
% Win = hamming(2*WinSize+1);
% Win = ones(1,2*WinSize+1);
Window = zeros(1,XPts);
Window(L-WinSize:L+WinSize) = 1;
parfor j = 1:XPts
    ESig_K_Cut(:,:,j) = ESig_K(:,:,j).*Window(j);
end

proj1 = sum(sum(abs(ESig_K_Cut),2),1);
proj1 = reshape(proj1,1,[]);
figure(2);
plot(kAxis,proj1);
xlabel('k');
ylabel('Signal amplitude');
saveas(gcf, strcat(Folder, 'WaveVectCut'), 'emf');

% IFT to spatial domain
ESig_X_Cut = ifft(ESig_K_Cut,[],3);
ESig_t = ESig_X_Cut(:,:,1);
ESig_t = reshape(ESig_t,Size1,Size1);
% ESig_t = interp2(tAxis,tauAxis',ESig_t,tAxis,tAxis');

% TI-FWM signal
FWM = sum(abs(ESig_t).^2,2);
figure(5);
plot(tauAxis1,FWM);
xlabel('\tau (ps)');
ylabel('TI-FWM Intensity');
saveas(gcf, strcat(Folder, 'TI-FWM'), 'emf');

%% Zero pad in time
PadSize = Size - Size1;
Pad1 = zeros(Size1,PadSize);
Pad2 = zeros(PadSize,Size);
ESig_t = cat(1,cat(2,ESig_t,Pad1),Pad2);
%% Define plot ranges
TimeMax = 20;    % pico-seconds
EgyMax = 3;         % meV (one side)
EgyC = 0;
TimeRange = [0 TimeMax 0 TimeMax];
EgyPlotRange = [EgyC-EgyMax EgyC+EgyMax -(EgyC+EgyMax) -(EgyC-EgyMax)];
NCont = 50;

VMax = max(max(abs(ESig_t)));
figure(3)
set(gcf, 'Units', 'inch');
set(gcf, 'position', [ 0.5 1 14 5 ]);
subplot(131);
contourf(tauAxis,tauAxis',abs(ESig_t),linspace(0,VMax,NCont),'LineStyle','none');
xlabel('t (ps)');
ylabel('\tau (ps)');
title('Absolute value');
axis(TimeRange);
colorbar ('Location','NorthOutside');
subplot(132);
contourf(tauAxis,tauAxis',real(ESig_t),linspace(-VMax,VMax,NCont),'LineStyle','none');
xlabel('t (ps)');
ylabel('\tau (ps)');
title('Real part');
axis(TimeRange);
colorbar ('Location','NorthOutside');
subplot(133);
contourf(tauAxis,tauAxis',imag(ESig_t),linspace(-VMax,VMax,NCont),'LineStyle','none');
xlabel('t (ps)');
ylabel('\tau (ps)');
title('Imaginary part');
axis(TimeRange);
colorbar ('Location','NorthOutside');
saveas(gcf, strcat(Folder, 'TimeDomain2D'), 'emf');

% ESig_t = ESig_t';
% NFin = 512;
% tauMax1 = 1E-11;     % Maximum scan time for tau
% tauStep1 = tauMax1/NFin;
% tauAxis1 = 0:tauStep1:tauMax1-tauStep1;
% tAxis1 = tauAxis1';
% ESig_t1 = interp2(tauAxis1,tAxis,ESig_t,tauAxis1,tAxis1);
% Evaluate frequency-domain signal
Spec_2D = fftshift(fft(fftshift(fft(ESig_t,[],2),2),[],1),1);
Spec_2D = flipud(Spec_2D);
% Spec_2D = interp2(Egy_t1,Egy_tau1,Spec_2D1,Egy_t,Egy_tau);
% Spec_2D = abs(Spec_2D);
VMax1 = max(max(abs(Spec_2D)));
% Ft = FAxis - FAxis(Size/2+1);
% Egy_t = EgyAxis;
% Ftau = -(FAxis - FAxis(Size/2))';
% Egy_tau = flipud(EgyAxis');
figure(4);
set(gcf, 'Units', 'inch');
set(gcf, 'position', [ 0.5 1 14 5 ]);
subplot(131);
contourf(Egy_t,Egy_tau,abs(Spec_2D),linspace(0,VMax1,NCont),'LineStyle','none');
line([EgyPlotRange(1) EgyPlotRange(2)], [EgyPlotRange(4) EgyPlotRange(3)], 'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
xlabel('\omega_t (meV)');
ylabel('\omega_\tau (meV)');
title('Absolute value');
axis(EgyPlotRange);
colorbar ('Location','NorthOutside');
subplot(132);
contourf(Egy_t,Egy_tau,real(Spec_2D),linspace(-VMax1,VMax1,NCont),'LineStyle','none');
line([EgyPlotRange(1) EgyPlotRange(2)], [EgyPlotRange(4) EgyPlotRange(3)], 'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
xlabel('\omega_t (meV)');
ylabel('\omega_\tau (meV)');
title('Real part');
axis(EgyPlotRange);
colorbar ('Location','NorthOutside');
subplot(133);
contourf(Egy_t,Egy_tau,imag(Spec_2D),linspace(-VMax1,VMax1,NCont),'LineStyle','none');
line([EgyPlotRange(1) EgyPlotRange(2)], [EgyPlotRange(4) EgyPlotRange(3)], 'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
xlabel('\omega_t (meV)');
ylabel('\omega_\tau (meV)');
title('Imaginary part');
axis(EgyPlotRange);
colorbar ('Location','NorthOutside');
saveas(gcf, strcat(Folder, 'FrequencyDomain2D'), 'emf');
