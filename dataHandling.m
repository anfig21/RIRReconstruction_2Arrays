function Data = dataHandling(Data)
%Data = dataHandling(Data) Processes the corresponding data of the given
%loudspeaker for a 2-spherical-array configuration.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% PRE-PROCESSING
Data.Ref.pos = vertcat(Data.Line1.pos,Data.Line2.pos,Data.Line3.pos);
Data.Ref.h = horzcat(Data.Line1.h,Data.Line2.h,Data.Line3.h);
Data.Ref.n = horzcat(Data.Line1.n,Data.Line2.n,Data.Line3.n);

Data = rmfield(Data,{'Line1','Line2','Line3'});

% Recwin up to T = 1s
Data.Ref.h = Data.Ref.h(1:Data.Nsamples,:);
Data.Ref.n = Data.Ref.n(1:2*Data.Nsamples,:);
Data.SphL.h = Data.SphL.h(1:Data.Nsamples,:);
Data.SphL.n = Data.SphL.n(1:2*Data.Nsamples,:);
Data.SphL.p = Data.SphL.p(1:2*Data.Nsamples,:);
Data.SphR.h = Data.SphR.h(1:Data.Nsamples,:);
Data.SphR.n = Data.SphR.n(1:2*Data.Nsamples,:);
Data.SphR.p = Data.SphR.p(1:2*Data.Nsamples,:);

%% SPHERICAL ARRAYS (Inner sphere: 156 Microphones: samples 155-end)
M = 78;
[idxL,Data.SphL.fNyq] = nearestNeighboursSph(Data.SphL.pos(155:end,:),M,Data.c);
[idxR,Data.SphR.fNyq] = nearestNeighboursSph(Data.SphR.pos(155:end,:),M,Data.c);

% Left sphere
Data.SphL.pos = Data.SphL.pos(idxL,:);
Data.SphL.h = Data.SphL.h(:,idxL);
Data.SphL.p = Data.SphL.p(:,idxL);
Data.SphL.M = size(Data.SphL.pos,1);

% Right sphere
Data.SphR.pos = Data.SphR.pos(idxR,:);
Data.SphR.h = Data.SphR.h(:,idxR);
Data.SphR.p = Data.SphR.p(:,idxR);
Data.SphR.M = size(Data.SphR.pos,1);

%% LOW-PASS FILTER
% Filter design
Fc = min([Data.SphL.fNyq Data.SphR.fNyq]);
lpFilt = designfilt('lowpassfir','PassbandFrequency',Fc, ...
            'StopbandFrequency',Fc*1.3,'PassbandRipple',0.5, ...
            'StopbandAttenuation',65,'DesignMethod','kaiserwin',...
            'SampleRate',Data.Fs);

% Reference
Data.Ref.h = filtfilt(lpFilt,Data.Ref.h);
[Data.Ref.H,~] = fftUniBi(Data.Ref.h);

% Left sphere
Data.SphL.h = filtfilt(lpFilt,Data.SphL.h); [Data.SphL.H,~] = fftUniBi(Data.SphL.h);
Data.SphL.p = filtfilt(lpFilt,Data.SphL.p); [Data.SphL.P,~] = fftUniBi(Data.SphL.p);
Data.SphL.n = filtfilt(lpFilt,Data.SphL.n); [Data.SphL.N,~] = fftUniBi(Data.SphL.n);

% Right sphere
Data.SphR.h = filtfilt(lpFilt,Data.SphR.h); [Data.SphR.H,~] = fftUniBi(Data.SphR.h);
Data.SphR.p = filtfilt(lpFilt,Data.SphR.p); [Data.SphR.P,~] = fftUniBi(Data.SphR.p);
Data.SphR.n = filtfilt(lpFilt,Data.SphR.n); [Data.SphR.N,~] = fftUniBi(Data.SphR.n);

end

