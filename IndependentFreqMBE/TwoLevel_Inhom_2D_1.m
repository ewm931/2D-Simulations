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

% Inhomogeneity is added. Write OBEs for all states in inhomogeneous
% distribution. All the imhomogeneity parameters are added in the
% function. In this program only the number of frequency points used is
% needed. An odd number of points are used to ensure symmetry. For multiple
% resonances, same number of points should be used. Will be easier to
% wrtie OBEs.
% parpool(12);
clear all;% clc;

meV2Hz = 241.79895E9*2*pi;      % hbar

%% Solve OBE's during tau

% Define inhomogeneous width (Make sure same values are used in the
% fucntion)
g_Inhom = 0.4*meV2Hz;   % Inhomogeneity (sigma)
FPts = 25;
InhomRange = 8*g_Inhom;
dInhom = InhomRange/(FPts-1);
InhomEgy = -0.5*InhomRange : dInhom : 0.5*InhomRange;
% Set amplitude for each frequency component (as dipole factor)
Amp = exp(-0.5*(InhomEgy./g_Inhom).^2);
Amp = Amp./sum(Amp);

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

% Define tau steps
Size = 200;
tauStep = 1E-13;
tauAxis = 0:tauStep:(Size-1)*tauStep;
FStep = 1/tauAxis(end);
EgyStep = FStep*2*pi/meV2Hz;
EgyAxis = -EgyStep*(Size/2):EgyStep:EgyStep*(Size/2-1);
Egy_t = EgyAxis;
Egy_tau = flipud(EgyAxis');

% FRange = EgyRange*meV2Hz/(2*pi);% Total energy range
tauMax = max(tauAxis);     % Maximum scan time for tau

T = 2E-13;          % Delay T

% Full time axis parameters
TimeStart = -7E-13;

% Define matrix for final Y value
% Y_All_F = zeros(Size,Size,3,FPts,XPts);
Y_All = zeros(Size,Size,3,XPts);
y0_F1 = zeros(3,FPts);
y0_F = reshape(y0_F1,3*FPts,1);

YTemp = zeros(Size,3,XPts);

% Define frequency mask
FMask = ones(Size,3,FPts);
parfor j = 1:FPts
    FMask(:,:,j) = FMask(:,:,j).*Amp(j);
end

% Define tolerance for ode
options = odeset('RelTol',1e-6,'AbsTol',1e-12);
tic;

for j = 1:Size
    tau = tauAxis(j);
    Pos2 = Pos1 + tau;
    Pos3 = Pos2 + T;
    TimeMax = Pos3 + tauMax + 1E-12;
    PosMat = [Pos1;Pos2;Pos3];
    
    parfor k = 1:XPts
        Phase = 2*pi*XAxis(k).*k_pulse;
        Phase(1) = Phase(1) + pi;
        XAmp = AmpX(k);
        y0 = cat(1,y0_F,Phase,PosMat,XAmp);
        [Time1,Y1] = ode23(@timeEvol_Inhom_2Lvl,[TimeStart,TimeMax],y0,options);
        Y1 = Y1(:,1:FPts*3);
        
        % Take points after the 3rd pulse is incident
        Time1 = Time1 - Pos3;
        Pos3_Idx = find(Time1 >= 0);
        Y1 = Y1(Pos3_Idx(1)-1:end,:);
        Time1 = Time1(Pos3_Idx(1)-1:end);
        
        Y = interp1(Time1',Y1,tauAxis');
        Y = reshape(Y,Size,3,FPts);
        
        Y_FMask = Y.*FMask;
        Y_Inhom = sum(Y_FMask,3);
        
        YTemp(:,:,k) = Y_Inhom;
    end
%     Y_All1 = zeros(Size,3,XPts);
%     parfor l = 1:FPts
%         Y_temp = Amp(l).*reshape(YTemp(:,:,l,:),Size,3,XPts);
%         Y_All1 = Y_All1 + Y_temp;
%     end
%     Y_All1 = reshape(Y_All1,1,Size,3,XPts);
    Y_All(j,:,:,:) = YTemp;
%     j
%     toc;
end
toc;
%%
% Save 3rd order polarization (sigma_10)
Pol3_X = Y_All(:,:,2,:) - 1i*Y_All(:,:,3,:);

ESig_X = 1i*Pol3_X;
dlmwrite('ESig_X1.dat',ESig_X,'\t');
toc;
%% Comment out from here for phase cycling
% ESig_K = fft(ESig_X,[],4);
% proj1 = sum(sum(abs(ESig_K),2),1);
% proj1 = reshape(proj1,1,[]);
% figure(1);
% plot(kAxis,proj1);
% dlmwrite('kProj.dat',proj1,'\t');
% dlmwrite('kAxis.dat',kAxis,'\t');
% 
% %%
% % Select out the 3rd order signal
% kSig = -k_pulse(1)+k_pulse(2)+k_pulse(3);
% kDiff = abs(kAxis - kSig);
% L = find(kDiff == min(kDiff));
% 
% % % Cut relevant points in k-space
% % NkPts = 2;
% ESig_K_Cut = zeros(size(ESig_K));
% % ESig_K_Cut(:,:,:,L-NkPts:L+NkPts) = ESig_K(:,:,:,L-NkPts:L+NkPts);
% 
% % Apply window to smooth the function
% % Win = hamming(5);
% Window = zeros(1,XPts);
% Window(L) = 1;
% parfor j = 1:XPts
%     ESig_K_Cut(:,:,:,j) = ESig_K(:,:,:,j).*Window(j);
% end
% 
% proj1 = sum(sum(abs(ESig_K_Cut),2),1);
% proj1 = reshape(proj1,1,[]);
% figure(6);
% plot(kAxis,proj1);
% 
% 
% % IFT to spatial domain
% ESig_X_Cut = ifft(ESig_K_Cut,[],4);
% ESig_t = ESig_X_Cut(:,:,:,1);
% ESig_t = reshape(ESig_t,Size,Size);
% dlmwrite('ESig_t.dat',ESig_t,'\t');
% 
% VMax = max(max(abs(ESig_t)));
% figure(5);
% contourf(tauAxis,tauAxis',abs(ESig_t'), linspace(0,VMax,50), 'LineStyle', 'none');
% axis square;
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
% 
% dlmwrite('Spec2D.dat',Spec_2D,'\t');
% 
% % % Close the parallel pool
% % poolobj = gcp('nocreate');
% % delete(poolobj);
% % figure(7);
% % contourf(real(Spec_2D),linspace(0,VMax1,50),'LineStyle','none');
