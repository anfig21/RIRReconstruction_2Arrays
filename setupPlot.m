function setupPlot(Data,flagS,flagH,flagT)
%setupPlot(Data,flagS,flagH,flagT) Plots the measurement setup, the
%frequency response and the RIR at the reference microphones.
%   Input:
%       - Data      : data structure.
%       - flagS     : plots setup
%                       'true'
%                       'false' (Default value)
%       - flagH     : plots frequency response
%                       'true'
%                       'false' (Default value)
%       - flagT     : plots RIR at reference line
%                       'true'
%                       'false' (Default value)
%
% Author: Antonio Figueroa Durán
% Date: November 2022

%% ERROR HANDLING
% plotFlag default value
if nargin < 1, error('setupPlot Error: Not enough input parameters.'), end
if nargin < 4, flagT = false; end
if nargin < 3, flagH = false; end
if nargin < 2, flagS = false; end

%% MAIN CODE
% Time vector
T = [25 45]*1e-3;       % Source near field
% T = [15 35]*1e-3;    % Source far field
t = Data.t(Data.t >= T(1) & Data.t <= T(2));

%% PLOT: SETUP
if flagS
    figure
    scatter3(Data.Source.pos(1),Data.Source.pos(2),Data.Source.pos(3),200,'filled'), hold on
    scatter3(Data.Ref.pos(:,1),Data.Ref.pos(:,2),Data.Ref.pos(:,3))
    scatter3(Data.SphL.pos(:,1),Data.SphL.pos(:,2),Data.SphL.pos(:,3))
    scatter3(Data.SphR.pos(:,1),Data.SphR.pos(:,2),Data.SphR.pos(:,3))
    drawRoom(Data.D(1),Data.D(2),Data.D(3)), axis equal
    axis([0 Data.D(1) 0 Data.D(2) 0 Data.D(3)])
    xlabel('$x$/m'), ylabel('$y$/m'), zlabel('$z$/m')
    legend('Source','Reference Line','Spherical Array 1','Spherical Array 2')
    applyAxisProperties(gca)
    applyLegendProperties(gcf)
end

%% PLOT: FREQUENCY RESPONSE
if flagH, plotFreqResponse(Data), end

%% PLOT: REFERENCE RIR
if flagT
    figure, hold on
    s = pcolor(Data.Ref.pos(:,1),t*1e3,Data.Ref.h(ismember(Data.t,t),:));
    set(s,'edgecolor','none')
    xlabel('$x$/m'), ylabel('$t$/ms')
    colormap hot
    c = colorbar;
%     caxis([-0.04 0.04])
    xline(min(Data.SphL.pos(:,1))), xline(max(Data.SphL.pos(:,1)))
    xline(min(Data.SphR.pos(:,1))), xline(max(Data.SphR.pos(:,1)))
    applyColorbarProperties(c,'$h(t)$/PaV$^{-1}$')
    applyAxisProperties(gca)
end

end

