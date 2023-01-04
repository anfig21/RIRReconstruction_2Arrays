%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  ROOM IMPULSE RESPONSE RECONSTRUCTION
%
% -------------------------------------------------------------------------
% Figueroa-Duran, Fernandez-Grande, "Reconstruction of room impulse
% responses over an extended spatial domain using block-sparse and kernel
% regression methods", International Congress of Acoustics, 2022
% -------------------------------------------------------------------------
%
% MOD: 2 spherical microphone arrays. Range estimated using DOA
%
% Antonio Figueroa-Duran
% anfig@dtu.dk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all

addpath(genpath(pwd))
addpath(genpath('../tools'))
% addpath(genpath('/Volumes/anfig/Data/room019/'))      % MacBook Pro
addpath(genpath('M:\Data'))                             % Windows 10

loadPlotParams

%% INITIAL PARAMETERS
Data.T = 1;                         % Pre-processing T = 1 s
Data.Fs = 48e3;
Data.t = 0:1/Data.Fs:Data.T-1/Data.Fs;
Data.D = [6.266,9.357,2.977];       % Room dimensions
Data.Nsamples = Data.Fs*Data.T;
Data.f = Data.Fs/Data.Nsamples*(0:Data.Nsamples/2-1);
Data.Temp = 23.2;
Data.Humi = 34.7;
Data.p0 = 100658;

[Data.rho,Data.c,~,Data.gamma,~,~,~,~,~,~] = amb2prop(Data.p0,Data.Temp,...
    Data.Humi,1);

Data.loudspeaker = 'rir_019_spk1.h5';
% Data.loudspeaker = 'rir_019_spk2.h5';

%% DATA ACQUISITION
Data = dataAcquisition(Data);

%% DATA HANDLING
Data = dataHandling(Data);

%% SETUP PLOT
% Flags
% - Setup
% - Frequency response
% - RIR reference line
setupPlot(Data,false,false,true);

%% ------------ DIRECT SOUND FIELD ------------ %%
% Windowing
Direct.T = 8.5*1e-3;      % Source near field
% Direct.T = 22*1e-3;     % Source far field

Direct = windowRIR(Data,0,Direct.T);

%% DOA Estimation
N = 1e3;                % Number of plane waves
DOAMethod = 'RLS';      % DOA Estimation Method
f = Data.f(1e2 <= Data.f & Data.f <= 1.5e3);
f = f(1:10:end);

Direct = earlyDOA(Data,Direct,N,f,DOAMethod);
clear N DOAMethod f

%% Range Estimation (DOA-based)
rMax = 4;

Direct.Range = earlyRange(Data,Direct,rMax);
clear rMax

%% Reconstruction
Rec.T = [5 7.5]*1e-3;
Rec.Direct = windowRIR(Data,Rec.T(1),Rec.T(2));

res = .5e-1;     % Spatial resolution (m)
L = 0.05;          % Cube side (m)
Rec.Direct.f = Data.f(0 <= Data.f & Data.f <= Data.SphL.fNyq);

Rec.Direct = clusterCoefficients(Data,Rec.Direct,Direct.Range,L,res);
clear res L

Rec.Direct = reconstructReflection(Data,Rec.Direct,Data.Ref.pos,true);

%% ------------ EARLY REFLECTIONS ------------ %%
% Peak detection using the frontmost microphone of each array
IR_L = Data.SphL.h(:,Data.SphL.pos(:,2) == max(Data.SphL.pos(:,2)));
IR_R = Data.SphR.h(:,Data.SphR.pos(:,2) == max(Data.SphR.pos(:,2)));

% Energy
ene_L = IR_L.^2;
ene_R = IR_R.^2;

% Peak detection
peaksL = findpeaksx(Data.t,ene_L,0,max(ene_L)*0.05,0,1,1);
peaksR = findpeaksx(Data.t,ene_R,0,max(ene_L)*0.05,0,1,1);

%% PLOT IR VS ENERGY
figure
subplot(221), plot(Data.t,IR_L), hold on
for ii = 1:size(peaksL,1), xline(peaksL(ii,2)), if(peaksL(ii,2)>25e-3), break,end, end, grid on
title('IR L'), xlim([5 25]*1e-3), applyAxisProperties(gca)

subplot(222), plot(Data.t,IR_R), hold on
for ii = 1:size(peaksR,1), xline(peaksR(ii,2)), if(peaksR(ii,2)>25e-3), break,end, end, grid on
title('IR R'), xlim([5 25]*1e-3), applyAxisProperties(gca)

subplot(223), plot(Data.t,ene_L), hold on
for ii = 1:size(peaksL,1), xline(peaksL(ii,2)), if(peaksL(ii,2)>25e-3), break,end, end, grid on
title('ene L'), xlim([5 25]*1e-3), applyAxisProperties(gca)

subplot(224), plot(Data.t,ene_R), hold on
for ii = 1:size(peaksR,1), xline(peaksR(ii,2)), if(peaksR(ii,2)>25e-3), break,end, end, grid on
title('ene R'), xlim([5 25]*1e-3), applyAxisProperties(gca)

%% INDEPENDANT DOA ESTIMATION
radius = 0.25;
margin = 2*radius/Data.c;     % Time to account for 



%% Windowing
Ref.T = [9 12]*1e-3;

Ref = windowRIR(Data,Ref.T(1),Ref.T(2));

%% DOA Estimation
N = 1e3;                % Number of plane waves
DOAMethod = 'RLS';      % DOA Estimation Method
f = Data.f(1 <= Data.f & Data.f <= 1.5e3);
f = f(1:10:end);

Ref = earlyDOA(Data,Ref,N,f,DOAMethod,true);
clear N DOAMethod f

%% Range Estimation (DOA-based)
rMax = 6;

Ref.Range = earlyRange(Data,Ref,rMax,true);
clear rMax

%% Reconstruction
Rec.T = [9 12]*1e-3;
Rec.Ref = windowRIR(Data,Rec.T(1),Rec.T(2));

res = .5e-1;     % Spatial resolution (m)
L = 0.1;          % Cube side (m)
Rec.Ref.f = Data.f(0 <= Data.f & Data.f <= Data.SphL.fNyq);

Rec.Ref = clusterCoefficients(Data,Rec.Ref,Ref.Range,L,res,true);
clear res L

%%
Rec.Ref = reconstructReflection(Data,Rec.Ref,Data.Ref.pos,true);

%% ------------ KERNEL RIDGE REGRESSION ------------ %%
Rec.T = [12 100]*1e-3;
Rec.KRR = windowRIR(Data,Rec.T(1),Rec.T(2));

Rec.KRR.t = Data.t(Data.t >= Rec.T(1) & Data.t <= Rec.T(2));
Rec.KRR.f = Data.f(0 <= Data.f & Data.f <= Data.SphL.fNyq);

Nf = numel(Rec.KRR.f);
Nr = size(Data.Ref.pos,1);
M = Data.SphL.M+Data.SphR.M;

% Euclidean distance matrices
d_MM = squareform(pdist([Data.SphL.pos; Data.SphR.pos]));
d_MNr = pdist2(Data.Ref.pos,[Data.SphL.pos; Data.SphR.pos])';

% Kernel Ridge Regression
k = 2*pi*Rec.KRR.f/Data.c;
Rec.KRR.H = zeros(Nr,numel(Data.f));

% sigma_2 = std(KRR.SphL.H(ismember(Data.f,KRR.f),:),[],2);

figure
for ii = 1:Nf
    p = [Data.SphL.H(Data.f==Rec.KRR.f(ii),:) Data.SphR.H(Data.f==Rec.KRR.f(ii),:)].';
    
    % Spherical Kernel
    K_MM = sinc(k(ii)*d_MM/pi);
    K_MNr = sinc(k(ii)*d_MNr/pi);
    
    % Cylindrical Kernel - Performs worse than Spherical Kernel
%     K_MM = sin(k(ii)*d_MM)./sqrt(k(ii)*d_MM);
%     K_MNr = sin(k(ii)*d_MNr)./sqrt(k(ii)*d_MNr);
%     
%     % Replace NaN or Inf
%     K_MM(k(ii)*d_MM==0) = 1;
%     K_MNr(k(ii)*d_MNr==0) = 1;
    
    % Normalisation
    K_MM = K_MM./vecnorm(K_MM);
    K_MNr = K_MNr./vecnorm(K_MNr);

    % Regularisation parameter selection method
    [U,s,V] = csvd(K_MM);
    %[lambda,~,~,~] = l_curve(U,s,p,'Tikh',[],[],false);
    [lambda,~,~] = quasiopt(U,s,p,'Tikh');
    
    Rec.KRR.H(:,Data.f==Rec.KRR.f(ii)) = K_MNr.'*((K_MM+lambda*eye(M))\p);
end

% Double-sided spectrum
KRR.P2 = [real(Rec.KRR.H(:,1)) Rec.KRR.H(:,2:end)/2];
KRR.P2 = [KRR.P2 flip(conj(KRR.P2),2)];
KRR.h = ifft(KRR.P2*Data.Nsamples,[],2,'symmetric').';

%%
figure, hold on
s = pcolor(Data.Ref.pos(:,1),Rec.KRR.t*1e3,KRR.h(ismember(Data.t,Rec.KRR.t),:));
set(s,'edgecolor','none')
xlabel('$x$/m'), ylabel('$t$/ms')
ylim(Rec.T*1e3)
colormap hot
c = colorbar;
% caxis([-0.04 0.04])
xline(min(Data.SphL.pos(:,1))), xline(max(Data.SphL.pos(:,1)))
xline(min(Data.SphR.pos(:,1))), xline(max(Data.SphR.pos(:,1)))
applyColorbarProperties(c,'$h(t)$/PaV$^{-1}$')
applyAxisProperties(gca)

%% ERROR
error2D = abs(KRR.h-Data.Ref.h)/norm(Data.Ref.h);
normError_Sph = norm(error2D(ismember(Data.t,Rec.KRR.t),:));

figure, hold on
s = pcolor(Data.Ref.pos(:,1),Rec.KRR.t*1e3,error2D(ismember(Data.t,Rec.KRR.t),:));
set(s,'edgecolor','none')
xlabel('$x$/m'), ylabel('$t$/ms')
ylim(Rec.T*1e3)
colormap hot
c = colorbar;
% caxis([-0.04 0.04])
xline(min(Data.SphL.pos(:,1))), xline(max(Data.SphL.pos(:,1)))
xline(min(Data.SphR.pos(:,1))), xline(max(Data.SphR.pos(:,1)))
applyColorbarProperties(c,'Relative Error')
applyAxisProperties(gca)
