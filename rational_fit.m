%Author : Yevgeniy Simonov
%Date : June 3, 2020
%Curtin University, 2020

close all 
clear all

Nexp = 3;
Npoles = 200;
freq_dev = 0.5e9; %how much you want to cutoff from start and end frequency

experiment = ["dataM/LOS.s2p", ... %EXP(1) Line-of-sight,3m
              "dataM/NOLOS.s2p", ... %EXP(2) No Line-of-sight, 3m
              "dataM/1R.s2p", ... %EXP(3) 1 Side Reflector, 3m
              "dataM/2R.s2p", ... %EXP(4) 2 Side Reflectors, 3m
              "dataM/3R.s2p", ... %EXP(5) 3 Side Reflectros, 3m
              "dataM/3R_NOLOS.s2p", ... %EXP(6) 3 Side Reflectors No Line-of-sight, 3m
              "dataM/3R_2D.s2p", ... %EXP(7) 3 Side Reflectors, 2 Diffraction Plates, 3m
              "dataM/3R_4D.s2p", ... %EXP(8) 3 Side Reflectors, 4 Diffraction Plates, 3m
              "dataM/4R_4D.s2p"]; %EXP(9) 4 Side Reflectors, 4 Diffraction Plates

S = sparameters(experiment(Nexp)).Parameters;
freq1 = sparameters(experiment(Nexp)).Frequencies;

S21 = squeeze(S(2,1,:));
TF1 = S21/2;

%remove lower frequency components that cause phase deviation
%the upper frequencies have to be removed too to keep  
%center frequency at the same position
freq_cutoff1 = freq1(1) + freq_dev;
freq_cutoff2 = freq1(end) - freq_dev;

%find the closest position to cutoff points
pos1 = (abs(freq1 - freq_cutoff1) < 0.01e9);
pos2 = (abs(freq1 - freq_cutoff2) < 0.01e9);
    
f1 = min(find(pos1 == 1));
f2 = max(find(pos2 == 1));

%create a new frequency array and update boundaries of transfer function
freq = freq1(f1:f2);
TF = TF1(f1:f2);

%fit transfer function into rational fit
fit = rationalfit(freq,TF,'NPoles',Npoles);
TF_fitted = freqresp(fit,freq);

%%
fCenter = 5e9;
period = 1/fCenter;
sampleTime = period/16;
signalLen = 8192;
ts = (0:signalLen-1)'*sampleTime; % 256 periods
ts_ns = ts * 1e9;

Ns = size(ts,1);
phi = unifrnd(-pi,pi,Ns,1);
input_AM = (1+unifrnd(-1,1,Ns,1)).*cos(2*pi*fCenter*ts);
input_PM = cos(phi).*cos(2*pi*fCenter*ts)-sin(phi).*sin(2*pi*fCenter*ts);
input_TH = randn(size(ts))+cos(2*pi*fCenter*ts);

output_AM = timeresp(fit,input_AM,sampleTime);
output_PM = timeresp(fit,input_PM,sampleTime);
output_TH = timeresp(fit,input_TH,sampleTime);

Npts = signalLen/8;

power_corr_AM_Tx = power_correlation(0.5*abs(input_AM).^2,Npts);
power_corr_PM_Tx = power_correlation(0.5*abs(input_PM).^2,Npts);
power_corr_TH_Tx = power_correlation(0.5*abs(input_TH).^2,Npts);

power_corr_AM_Rx = power_correlation(0.5*abs(output_AM).^2,Npts);
power_corr_PM_Rx = power_correlation(0.5*abs(output_PM).^2,Npts);
power_corr_TH_Rx = power_correlation(0.5*abs(output_TH).^2,Npts);

[resp,t]=impulse(fit,1e-11,5e4);

%% Channel Transfer Function Magnitude and Phase
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
semilogy(freq/1e9,abs(TF),freq/1e9,abs(TF_fitted),'--','LineWidth',2)
xlabel('Frequency (GHz)')
ylabel('Magnitude')
legend('data','fit')
title('Channel Transfer Function and the Rational Fit')
grid on

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
plot(freq/1e9,180/pi*atan2(imag(TF),real(TF)))
grid on
xlabel('Frequency, GHz')
ylabel('Phase (Deg)')
title('Phase of the Transfer Function (Data)')
subplot(2,1,2)
plot(freq/1e9,unwrap(angle(TF)),...
     freq/1e9,unwrap(angle(TF_fitted)), '--','LineWidth',2)
grid on
xlabel('Frequency, GHz')
ylabel('Unwrapper Phase (Deg)')
legend('data','fit')

%% Channel Impulse Response
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(t*1e9,resp)
grid on
title('Impulse Response')
xlabel('Time,ns')
ylabel('Magnitude, linear')
xlim([0 25])

%% Signals at the transmitter
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(ts_ns,input_AM)
hold on
plot(ts_ns,input_PM)
hold on
plot(ts_ns,input_TH)
grid on
xlabel('Time,ns')
ylabel('Magnitude, V')
legend('AM','PM','Thermal')
xlim([0 2.5])
title('Transmitted Signals')

%% Signals at the receiver
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(ts_ns,output_AM)
hold on
plot(ts_ns,output_PM)
hold on
plot(ts_ns,output_TH)
grid on
xlabel('Time,ns')
ylabel('Magnitude, V')
legend('AM','PM','Thermal')
xlim([10 16])
title('Signals at the Receiver')

%% Power Correlation at Tx
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(ts(1:Npts)*1e9,power_corr_AM_Tx)
hold on
plot(ts(1:Npts)*1e9,power_corr_PM_Tx)
hold on
plot(ts(1:Npts)*1e9,power_corr_TH_Tx)
grid on
xlim([0 3])
legend('AM','PM','Thermal')
title('Power Correlation at the Transmitter')
xlabel('Time Delay (ns)')

%% Power Correlation at Rx
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(ts(1:Npts)*1e9,power_corr_AM_Rx)
hold on
plot(ts(1:Npts)*1e9,power_corr_PM_Rx)
hold on
plot(ts(1:Npts)*1e9,power_corr_TH_Rx)
grid on
xlim([0 3])
legend('AM','PM','Thermal')
title('Power Correlation at the Receiver')
xlabel('Time Delay (ns)')

%% Frequency Response
NFFT = 2^nextpow2(signalLen);
samplingFreq = 1/sampleTime;

df = freq(2)-freq(1);
PSD_input_AM = abs(fft(input_AM,NFFT)).^2/signalLen;
PSD_input_PM = abs(fft(input_PM,NFFT)).^2/signalLen;
PSD_input_TH = abs(fft(input_TH,NFFT)).^2/signalLen;

PSD_output_AM = abs(fft(output_AM,NFFT)).^2/signalLen;
PSD_output_PM = abs(fft(output_PM,NFFT)).^2/signalLen;
PSD_output_TH = abs(fft(output_TH,NFFT)).^2/signalLen;

f = samplingFreq/2*linspace(0,1,NFFT/2+1)';

%%
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(f/1e9,20*log10(PSD_input_AM(1:NFFT/2+1)),'-*','MarkerSize',3)
hold on
plot(f/1e9,20*log10(PSD_input_PM(1:NFFT/2+1)),'r','LineWidth',1)
hold on
plot(f/1e9,20*log10(PSD_input_TH(1:NFFT/2+1)),'-.')
legend('AM','PM','Thermal')
title('Power Spectral Density of the Transmitted Signals (Unnormalized)')
ylabel('Power / Frequency (dB/Hz)')
xlabel('Frequency, GHz')
grid on

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(f/1e9,20*log10(PSD_output_AM(1:NFFT/2+1)),'-*','MarkerSize',3)
hold on
plot(f/1e9,20*log10(PSD_output_PM(1:NFFT/2+1)),'r','LineWidth',1)
hold on
plot(f/1e9,20*log10(PSD_output_TH(1:NFFT/2+1)),'-.')
legend('AM','PM','Thermal')
title('Power Spectral Density of the Received Signals (Unnormalized)')
ylabel('Power / Frequency (dB/Hz)')
xlabel('Frequency, GHz')
grid on

function c = power_correlation(P,Npts)
    c = zeros(Npts,1);
    Pm = mean(P);
    Pmsqr = Pm^2;
    Np = size(P,1);
    for i=1:Npts
       Q = zeros(Np,1);
       Q(i+1:end)=P(1:end-i);
       Qm = mean(Q);
       c(i) = mean(P.*Q)/Qm/Pm; 
    end
end