function [] = plotFreqResponse(Data)
%plotFreqResponse(Data) Plots the measurement setup, the magnitude of the
%transfer function and the frequency response of signal and noise.
%   Input:
%       - Data      : data structure.
%
% Author: Antonio Figueroa Durán
% Date: November 2022

mic = 40;

% Frequency response
figure, plot(Data.f*1e-3,20*log10(abs([Data.SphL.H(:,mic) Data.SphR.H(:,mic)]))), grid on
xlabel('$f$/kHz'), ylabel('$|H(j\omega)|$/dB')
legend('Mic 40@Sph1','Mic 40@SphL')
applyAxisProperties(gca)
applyLegendProperties(gcf)

figure
subplot(221)
plot(Data.f*1e-3,20*log10(abs([Data.SphL.P(1:2:end,mic) Data.SphR.P(1:2:end,mic)])/20e-6)), grid on
xlabel('$f$/kHz'), ylabel('$|P(j\omega)|$/dB SPL')
legend('Mic 40@Sph1','Mic 40@SphL')
applyAxisProperties(gca)
applyLegendProperties(gcf)
ylim([-5 70])

subplot(222)
plot(Data.f,20*log10(abs([Data.SphL.P(1:2:end,mic) Data.SphR.P(1:2:end,mic)])/20e-6)), grid on
xlabel('$f$/Hz'), ylabel('$|P(j\omega)|$/dB SPL')
legend('Mic 40@Sph1','Mic 40@SphL')
applyAxisProperties(gca)
applyLegendProperties(gcf)
axis([0 4e3 -5 70])

subplot(223)
plot(Data.f*1e-3,20*log10(abs([Data.SphL.N(1:2:end,:) Data.SphR.N(1:2:end,:)])/20e-6)), grid on
xlabel('$f$/kHz'), ylabel('$|N(j\omega)|$/dB SPL')
applyAxisProperties(gca)
ylim([-5 70])

subplot(224)
plot(Data.f,20*log10(abs([Data.SphL.N(1:2:end,:) Data.SphR.N(1:2:end,:)])/20e-6)), grid on
xlabel('$f$/Hz'), ylabel('$|N(j\omega)|$/dB SPL')
applyAxisProperties(gca)
axis([0 4e3 -5 70])

end

