function Early = earlyDOA(Data,Early,N,f,DOAMethod,plotFlag)
%Direct = earlyDOA(Data,Early,N,f,DOAMethod,plotFlag)
%Estimates the Direction-of-Arrival (DOA) of the early part of the RIR
%using optimisation methods.
%   Input:
%       - Data      : raw data. Structure
%       - Early     : early part of the RIR (sparse). Structure
%       - N         : number of plane waves. Integer
%       - f         : frequency span. 1 x Nf
%       - DOAMethod : specifies the optimisation method:
%                       'RLS' - Tikhonov with L-curve
%                       'CS' - Compressive Sensing
%       - plotFlag  : 'true' to plot DOA estimation
%                     'false' (Default value)
%   Output:
%       - Early     : DOA Estimation. Structure
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
if nargin < 5, error('directSoundDOA Error: Not enough input parameters.'), end
if N <= 0, error('directSoundDOA Error: The number of plane waves must be a positive integer.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE
% Dictionary of plane waves
[HL,ukL] = dictionary(Data.c,f,Data.SphL.pos',N);
[HR,ukR] = dictionary(Data.c,f,Data.SphR.pos',N);
disp('Plane Wave Dictionary... OK')

% Optimisation problem
switch DOAMethod
    case 'RLS'
        % DOA Estimation via Regularised Least-Squares
        Early.DOA.SphL = dirDOA_RLS(Early.SphL.H(ismember(Data.f,f),:),HL,ukL);
        Early.DOA.SphR = dirDOA_RLS(Early.SphR.H(ismember(Data.f,f),:),HR,ukR);
    case 'CS'
        % DOA Estimation via Compressive Sensing
        Early.DOA.SphL = dirDOA_CS(Early.SphL.H(ismember(Data.f,f),:),...
            Early.SphL.Nnorm(ismember(Data.f,f)),HL,ukL);
        Early.DOA.SphR = dirDOA_CS(Early.SphR.H(ismember(Data.f,f),:),...
            Early.SphR.Nnorm(ismember(Data.f,f)),HR,ukR);
    otherwise
        error('earlyDOA Error: DOA estimation method not valid.')
end

% Include frequency vector in structure
Early.DOA.f = f;

% Mode
Early.DOA.ModeL = mode(Early.DOA.SphL,2);
Early.DOA.ModeR = mode(Early.DOA.SphR,2);

%% PLOT
if plotFlag    
    % 3-D Estimation
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on    % Ref line
    scatter3(Data.SphL.pos(:,1),Data.SphL.pos(:,2),Data.SphL.pos(:,3))          % Sph L
    scatter3(Data.SphR.pos(:,1),Data.SphR.pos(:,2),Data.SphR.pos(:,3))          % Sph R
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled')
    quiver3(Data.SphL.R0(1),Data.SphL.R0(2),Data.SphL.R0(3),...
        Early.DOA.ModeL(1),Early.DOA.ModeL(2),Early.DOA.ModeL(3),1.8,'Linewidth',4)
    quiver3(Data.SphR.R0(1),Data.SphR.R0(2),Data.SphR.R0(3),...
        Early.DOA.ModeR(1),Early.DOA.ModeR(2),Early.DOA.ModeR(3),1.8,'Linewidth',4)
    
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    xlabel('$x$/m'), ylabel('$y$/m'), zlabel('$z$/m')
    legend('Reference Line','Sphere@Left','Sphere@Right','Source',...
        'DOA Estimation@L','DOA Estimation@R')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

