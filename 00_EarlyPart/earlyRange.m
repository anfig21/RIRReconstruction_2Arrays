function range = earlyRange(Data,Early,rMax,plotFlag)
%Early = earlyRange(Data,Early,res,rMinMax,plotFlag)
%Estimates the Range of the apparent origin of the corresponding reflection
%based on DOA triangulation.
%   Input:
%       - Data          : raw data. Structure
%       - Early         : early part of the RIR. Structure
%       - rMax          : maximum distances from each mic array. Scalar
%       - plotFlag      : 'true' to plot Range estimation
%                         'false' (Default value)
%   Output:
%       - Early         : range estimation. 3 x 1
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
if nargin < 3, error('earlyRange Error: Not enough input parameters.'), end
if ~isfield(Early,'DOA'), error('earlyRange Error: DOA has not been estimated.'), end
if rMax < 0, error('earlyRange Error: rMax must be a positive scalar.'), end
if nargin < 4, plotFlag = false; end

%% MAIN CODE
kvect = -1:1e-3:rMax;      % Line parameter

yL = @(k) Early.DOA.ModeL*k+Data.SphL.R0';
yR = @(k) Early.DOA.ModeR*k+Data.SphR.R0';

yLvect = yL(kvect);
yRvect = yR(kvect);

% Direction of line connecting the closest points
n = cross(Early.DOA.ModeL,Early.DOA.ModeR);

% Parameter for closest points on each line
kL = dot(cross(Early.DOA.ModeR,n),(Data.SphR.R0-Data.SphL.R0))/dot(n,n);
kR = dot(cross(Early.DOA.ModeL,n),(Data.SphR.R0-Data.SphL.R0))/dot(n,n);

% Closest points on each line
yLClosest = yL(kL);
yRClosest = yR(kR);

% Source estimation
range = (yLClosest+yRClosest)/2;

%% PLOT
if plotFlag
    figure
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3)), hold on    % Ref line
    scatter3(Data.SphL.pos(:,1),Data.SphL.pos(:,2),Data.SphL.pos(:,3))          % Sph L
    scatter3(Data.SphR.pos(:,1),Data.SphR.pos(:,2),Data.SphR.pos(:,3))          % Sph R
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),150,'filled')

    plot3(yLvect(1,:),yLvect(2,:),yLvect(3,:))
    plot3(yRvect(1,:),yRvect(2,:),yRvect(3,:))

    scatter3(range(1),range(2),range(3),200,'filled')

    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    xlabel('x in m'), ylabel('y in m'), zlabel('z in m')
    legend('Reference Line','Sphere@Left','Sphere@Right','Source',...
        'DOA Estimation@L','DOA Estimation@R','Range Estimation')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

end

