clear;
clc;
close all;

% Time axis

t = linspace(0,150,500);

% Ideal magnetic anomaly

peakAmplitude = 1.4;
peakCentre = 80;
peakWidth = 10;

idealSignal = peakAmplitude .* ...
    exp(-(t-peakCentre).^2/(2*peakWidth^2));

figure('Color','w')

plot(t,idealSignal,...
    'LineWidth',2)

xlabel('Time/s')
ylabel('Signal')

title('Magnetic Anomaly Signal in Ideal Condition')

xlim([0 150])
ylim([0 1.5])

grid on
box on

set(gca,...
    'FontSize',11,...
    'LineWidth',1.2)

% Actual measurement

baseline = 4.7572e4;

backgroundTrend = -2.1*(t/150);

slowVariation = ...
      0.22*sin(2*pi*t/170) ...
    + 0.12*cos(2*pi*t/80);

noise = 0.28*randn(size(t));

% Make anomaly much weaker than the ideal case

actualSignal = baseline ...
             + backgroundTrend ...
             + slowVariation ...
             + 1.25*idealSignal ...
             + noise;

figure('Color','w')

plot(t,actualSignal,...
    'LineWidth',1.5)

xlabel('Time/s')
ylabel('B/nT')

title('Magnetic Anomaly Signal in Actual Condition')

xlim([0 150])

grid on
box on

set(gca,...
    'FontSize',11,...
    'LineWidth',1.2)