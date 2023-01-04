function Rec = clusterCoefficients(Data,Rec,R0,L,res,plotFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% ERROR HANDLING
if nargin < 6, plotFlag = false;    % plotFlag default value
elseif nargin < 5, error('clusterCoefficients Error: Not enough input parameters.'), end

%% MAIN CODE
% Combine data from two spheres
pos = [Data.SphL.pos; Data.SphR.pos];
P = [Rec.SphL.H(ismember(Data.f,Rec.f),:) Rec.SphR.H(ismember(Data.f,Rec.f),:)].';
Nnorm = mean([Rec.SphL.Nnorm Rec.SphR.Nnorm],2);

% Point sources - cubic grid
l = -L/2:res:L/2;
[X,Y,Z] = meshgrid(l,l,l);
Rec.rs = [X(:) Y(:) Z(:)]+R0(:)';

% figure, scatter3(rs(:,1),rs(:,2),rs(:,3)), axis equal

N = size(Rec.rs,1);             % Number of point sources
M = size(pos,1);                % Number of microphones
Nf = length(Rec.f);             % Number of frequency bins

k = (2*pi*Rec.f)/Data.c;        % Propagation vector

% Euclidean distance matrix
d = nan(M,N);                   % mics x point sources
for mm = 1:M
    d(mm,:) = vecnorm(repmat(pos(mm,:),N,1)-Rec.rs,2,2);
end

% Coefficient estimation
NoiseMargin = 10;
x = zeros(N,Nf);

c = parcluster('local');
parpool(c,2);
parfor ii = 1:Nf
    H = (1./d).*exp(-1i*d*k(ii));       % Dictionary (point sources)
    
    % CS solution
    epsilon = 10^(NoiseMargin/20)*Nnorm(ii);            
    x(:,ii) = compressiveSensing(H,P(:,ii),epsilon);
            
    disp(strcat("Reconstructing reflection... ",string(ii-1),'/',string(Nf-1)," Hz"))  
end

Rec.x = x;
delete(gcp('nocreate'))

%% PLOT
% Coefficients
if plotFlag
    figure
    plot(Rec.f,abs(Rec.x).'), grid on
    xlabel('$f$/Hz'), ylabel('$|\hat{x}|$')
    applyAxisProperties(gca)
end

end

