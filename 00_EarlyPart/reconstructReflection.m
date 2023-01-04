function Rec = reconstructReflection(Data,Rec,r,plotFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% ERROR HANDLING
if nargin < 4, plotFlag = false;    % plotFlag default value
elseif nargin < 3, error('reconstructReflection Error: Not enough input parameters.'), end

%% RECONSTRUCTION
R = size(r,1);
N = size(Rec.x,1);
Nf = length(Rec.f);

k = (2*pi*Rec.f)/Data.c;        % Propagation vector

d = nan(N,R);                   % mics x point sources
for nn = 1:N
    d(nn,:) = vecnorm(repmat(Rec.rs(nn,:),R,1)-r,2,2);
end

% Dictionary
Rec.P = zeros(R,length(Data.f));
for ii = 1:Nf
    %H = (1./d).*exp(-1i*d*k(ii));       % Dictionary (point sources)
    H = exp(-1i*d*k(ii));
    Rec.P(:,Data.f==Rec.f(ii)) = H.'*Rec.x(:,ii);
end

% Double-sided spectrum
P2 = [real(Rec.P(:,1)) Rec.P(:,2:end)/2];
P2 = [P2 flip(conj(P2),2)];
Rec.h = ifft(P2*Data.Nsamples,[],2,'symmetric').';

%% PLOT: REFERENCE RIR
if plotFlag
    T = [5 25]*1e-3;       % Source near field
    t = Data.t(Data.t >= T(1) & Data.t <= T(2));

    figure, hold on
    s = pcolor(r(:,1),t*1e3,Rec.h(ismember(Data.t,t),:));
    set(s,'edgecolor','none')
    xlabel('$x$/m'), ylabel('$t$/ms')
    colormap hot
    c = colorbar;
    caxis([-0.04 0.04])
    xline(min(Data.SphL.pos(:,1))), xline(max(Data.SphL.pos(:,1)))
    xline(min(Data.SphR.pos(:,1))), xline(max(Data.SphR.pos(:,1)))
    applyColorbarProperties(c,'$h(t)$/PaV$^{-1}$')
    applyAxisProperties(gca)
end

end

