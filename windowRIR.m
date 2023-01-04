function Structure = windowRIR(Data,Tini,Tfin,plotFlag)
%Structure = windowRIR(Data,Tini,Tfin,plotFlag) Applies Hanning window to
%RIR between time given by Tini and Tfin. Obtains the frequency response
%and the noise spectrum. Estimates the noise power after windowing.
%   Input:
%       - Data      : raw data. Structure
%       - Tini      : initial time. Scalar
%       - Tfin      : end time. Scalar
%       - plotFlag  : 'true'  - RIR in time domain
%                     'false' - Default value
%   Output:
%       - Structure : structure with windowed data
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 4, plotFlag = false;
elseif nargin < 3, error('windowRIR Error: Not enough input parameters.'), end

%% MAIN CODE
tHann = 1e-3;
NHann = Data.Fs*tHann/2;
h = hann(2*NHann);
hHalf1 = h(1:end/2);
hHalf2 = h(end/2+1:end);

Structure.T = [Tini Tfin];
Structure.N = floor(Data.Fs*Structure.T);
Structure.Nsamples = Structure.N(2)-Structure.N(1);

% Windowing - 1/2 Hanning (hann) window on each side
w = vertcat(zeros(Structure.N(1),Data.SphL.M),...
    repmat(hHalf1,1,Data.SphL.M),...
    ones(Structure.Nsamples-2*NHann,Data.SphL.M),...
    repmat(hHalf2,1,Data.SphL.M),...
    zeros(Data.Nsamples-Structure.Nsamples-Structure.N(1),Data.SphL.M));

% Windowing Sphere Left
Structure.SphL.h = w.*Data.SphL.h;
Structure.SphL.n = w(:,size(Data.SphL.n,2)).*Data.SphL.n(1:2:end,:);
[Structure.SphL.H,~] = fftUniBi(Structure.SphL.h);
[Structure.SphL.N,~] = fftUniBi(Structure.SphL.n);

% Windowing Sphere Right
Structure.SphR.h = w.*Data.SphR.h;
Structure.SphR.n = w(:,size(Data.SphR.n,2)).*Data.SphR.n(1:2:end,:);
[Structure.SphR.H,~] = fftUniBi(Structure.SphR.h);
[Structure.SphR.N,~] = fftUniBi(Structure.SphR.n);

% Noise norm
Structure.SphL.Nnorm = mean(abs(Structure.SphL.N),2);
Structure.SphR.Nnorm = mean(abs(Structure.SphR.N),2);

%% PLOT RIR
if plotFlag
    figure
    % Sphere Left
    subplot(221), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Data.SphL.h), grid on
    xlim([Tini*1e3-1 Tfin*1e3+1]), title('Sphere@Left'), ylabel('$h(t)$')
    applyAxisProperties(gca)
    subplot(223), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Structure.SphL.h), grid on
    xlim([Tini*1e3-1 Tfin*1e3+1]), ylabel('$h_w(t)$'), xlabel('$t$/ms')
    applyAxisProperties(gca)
    
    subplot(222), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Data.SphR.h), grid on
    xlim([Tini*1e3-1 Tfin*1e3+1]), title('Sphere@Right'), ylabel('$h(t)$')
    applyAxisProperties(gca)
    subplot(224), plot((0:Data.Nsamples-1)*1e3/Data.Fs,Structure.SphR.h), grid on
    xlim([Tini*1e3-1 Tfin*1e3+1]), ylabel('$h_w(t)$'), xlabel('$t$/ms')
    applyAxisProperties(gca)
end

end

