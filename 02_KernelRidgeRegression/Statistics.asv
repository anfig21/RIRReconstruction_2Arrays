clc, clear, close all
%% Mean square pressure
N = 5e3;    % Number of waves
M = 1e3;    % Number of fields
rng(0), phi = rand(N,M)*2*pi;

p_mean = abs(sum(exp(1j*phi),1)).^2/N;
p_amp = sqrt(p_mean);
SPL = 20*log10(p_amp/20e-6);

disp(strcat("Mean: ",string(mean(p_mean))))
disp(strcat("STD: ",string(std(p_mean))))
disp(strcat("Relative STD amplitude: ",string(std(p_amp)/mean(p_amp))))
disp(strcat("STD SPL: ",string(std(SPL))))

%% Spatial Correlation
O = 2e2;    % Number of points in the radial direction
N = 3e2;    % Number of waves
M = 1e3;    % Number of fields
rng(0)
phi = rand(N,M)*2*pi;
kr = linspace(0,5*pi,O);
theta = asin(2*(rand(N,M)-0.5))+pi/2;

p_xcorr = nan(O,M);
for mm = 1:M
    p_xcorr(:,mm) = sum(exp(1j*phi(:,mm)))*conj(sum(exp(1j*(phi(:,mm)-kr.*cos(theta(:,mm))))));
end
p_xcorr = p_xcorr/N;
p_mean = abs(sum(exp(1j*phi),1)).^2/N;
p_nxcorr = mean(p_xcorr,2)./mean(p_mean);

%%
figure, plot(kr,real(p_nxcorr),kr,sinc(kr/pi)), legend('Stats','Sinc'), grid on

