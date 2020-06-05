%Author: Yevgeniy Simonov
%Date: June, 2020
%Simulator of Rayleigh Fading Channel
%Ref: Rappaport: Wireless Communications - Principles and Practice

clear all
close all

%======== DEFINITION OF VARIABLES =========
N = 2^14; %number of frequency domain points (has to be a power of 2)
fc = 5e9; %carrier frequency in Hz
vmob = 0.01; %m/s speed of the mobile
%==========================================

Np = N/2; %number of one-sided frequency points
wc = 2*pi*fc; 
c = 3e8; %speed of RF wave in the air in m/s
fDmax = vmob*fc/c; %maximum Doppler shift in Hz
df = 2*fDmax/(N-1); %frequency resolution
T = 1/df; %duration of fading waveform
f = linspace(fc-fDmax, fc+fDmax, N).'; %array of frequencies
t = linspace(0, T, N).'; %time array

%positive frequency components 
Xp1 = randn(Np,1) + 1i * randn(Np,1);
Xp2 = randn(Np,1) - 1i * randn(Np,1);
    
%negative frequency components
Xn1 = conj(flipud(Xp1));
Xn2 = conj(flipud(Xp2));
    
%Frequency-domain in-phase and quadrature 
%components of the lowpass signal
I = vertcat(Xn1, Xp1); %in-phase 
Q = vertcat(Xn2, Xp2); %quadrature
    
SEz = zeros(N, 1); %Doppler spectral density
SEz(2:end-1) = 1.5./(pi*fDmax*sqrt(1-((f(2:end-1)-fc)/fDmax).^2));
SEz(1) = 2*SEz(2)-SEz(3);
SEz(end) = 2*SEz(end-1)-SEz(end-2);

%multiply sqrt of Doppler Spectrum with in-phase and qudrature
IS = I.*sqrt(SEz);
QS = Q.*sqrt(SEz);
    
%compute inverse fourier transform of the frequency-domain signals
It = ifft(IS);
Qt = ifft(QS);

Envelope = sqrt(abs(It).^2 + abs(Qt).^2);
Phase = atan2(abs(Qt),abs(It))*180/pi;

%%
gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);

subplot(2,1,1)
plot(t,20*log10(Envelope))
grid on
xlabel('Time, s')
ylabel('Envelope Amplitude, dB')
title('Characteristics of Rayleigh Fading Channel')
xlim([min(t) max(t)])

subplot(2,1,2)
plot(t,unwrap(Phase))
grid on
xlabel('Time, s')
ylabel('Unwrapped Phase, Deg')
xlim([min(t) max(t)])

%% Amplitude Distribution
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
p = histogram(Envelope,'Normalization','probability','BinWidth',0.001); 
title('Probability Density Function of Rayleigh Channel')
grid on
xlabel('Envelope Magnitude')