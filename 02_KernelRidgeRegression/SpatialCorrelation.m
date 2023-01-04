%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Advanced Acoustics
% Exercise 3: The diffuse sound field
%   Spatial correlation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all

%% INITIAL PARAMETERS
N = 100;            % Number of incident waves
numAv = 1E3;        % Number of repetitions
R = 2E2;            % Spatial resolution

kr = linspace(0,5*pi,R);

%% MONTE CARLO SIMULATION
p1p2 = zeros(R,1);
pRMS = zeros(R,1);

% Phi and Theta matrix: numAv x N
phi = rand(numAv,N)*2*pi;
theta = asin(2*(rand(numAv,N)-0.5))+pi/2;

% Spatial correlation
c = waitbar(0,'Loading...0%','Name','Calculation running');
for ii = 1:R
    p1p2(ii) = mean(sum(exp(1i.*phi),2).*conj(sum(exp(1i*(phi-kr(ii)*cos(theta))),2)));
    waitbar(ii/R,c,sprintf('Loading... %.f%%',100*ii/R));
end
delete(c)

% Mean square pressure
pRMS(:) = mean(abs(sum(exp(1i.*phi),2)).^2);


%% REFERENCE
ref = sin(kr)./kr;
ref(ref == 0) = 1;

%% PLOT DATA
figure
plot(kr,ref,'--','Linewidth',2), hold on, grid on
plot(kr,real(p1p2)./pRMS,'Linewidth',2)

set(get(gca,'XAxis'),'Fontsize',23)
set(get(gca,'YAxis'),'Fontsize',23)
xticks(0:pi:5*pi), xlim([0 5*pi])
xticklabels({'0','\pi','2\pi','3\pi','4\pi','5\pi'})

xlabel('\textbf{kr}','Interpreter','Latex','Fontsize',28)
ylabel('\textbf{Correlation Coefficient}','Interpreter','Latex','Fontsize',28)
legend('\textbf{Reference $sin(kr)/kr$}',...
    '\textbf{E$\{\hat{p}_1\hat{p}_2^*\}/$E$\{|\hat{p}|^2\}$}',...
    'Location','Best','Interpreter','Latex','Fontsize',28)

