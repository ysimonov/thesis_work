clear all
close all

Nt = 3000;
Nf = 12000;
time_multiplier = 80;
Nrays = 10000;
include_LOS = false;

Ro = 50;
f1C = 0.1e9; %lower 3dB cutoff
f2C = 15e9; %upper 3dB cutoff

%Butterworth bandpass filter
Ls = (Ro / (pi*(f2C - f1C)))/2;
Cs = 2*(f2C - f1C)/(4*pi*Ro*f2C*f1C);

Lp = 2*Ro*(f2C - f1C)/(4*pi*f2C*f1C);
Cp = (1/(pi*Ro*(f2C - f1C)))/2;

ckt = circuit('butterworthBPF'); 

add(ckt,[3 2],inductor(Ls))
add(ckt,[4 3],capacitor(Cs))
add(ckt,[5 4],capacitor(Cs))
add(ckt,[6 5],inductor(Ls))

add(ckt,[4 1],capacitor(Cp))
add(ckt,[4 1],inductor(Lp))
add(ckt,[4 1],inductor(Lp))
add(ckt,[4 1],capacitor(Cp))

freq = linspace(1e9,9e9,Nf).';
omega = 2*pi*freq;

setports(ckt,[2 1],[6 1])
S = sparameters(ckt,freq);

tfS = s2tf(S);
df = freq(2)-freq(1);
Tspan = 1/df;
time = linspace(-10e-9, 15e-9, Nt).';
dt = time(2)-time(1);
sample_rate_multiplier = Tspan / (time(end)-time(1)) * Nt / Nf;

CIR_LOS = IFT(tfS,omega,time);

positions = sort(randi(Nt,Nrays,1),'ascend');
attenuation_coefficients = rand(Nrays,1); %sort(rand(Nrays,1),'descend'); 

CIR_Nrays = zeros(Nt,1);
for i=1:Nrays
    dCIR = zeros(Nt,1);
    dCIR=circshift(CIR_LOS,positions(i));
    CIR_Nrays = CIR_Nrays + attenuation_coefficients(i)*dCIR;
end
if(include_LOS) 
    CIR_Nrays = CIR_Nrays + 100*CIR_LOS;
end
CIR_Nrays = CIR_Nrays / sum(CIR_Nrays); 

TF_LOS = FT(CIR_LOS,omega,time,sample_rate_multiplier);
TF_Nrays = FT(CIR_Nrays,omega,time,sample_rate_multiplier);

%% create test signals, pass through CIR and perform power correlation
Ns = time_multiplier * Nt;
ts = linspace(0,Ns,Ns)'*dt;

fc = 5e9;
wc = 2*pi*fc;

%narrowband phase modulation
phase = unwrap(wc*ts) + unwrap(2*pi*rand(Ns,1));
PM_input = real(exp(1i*phase));
AM_input = (1+unifrnd(-1,1,Ns,1)).*cos(2*pi*fc*ts);
TH_input = randn(size(ts))+cos(2*pi*fc*ts);

output_AM_LOS = convolution(AM_input,real(CIR_LOS));
output_PM_LOS = convolution(PM_input,real(CIR_LOS));
output_TH_LOS = convolution(TH_input,real(CIR_LOS));

output_AM_Nrays = convolution(AM_input,real(CIR_Nrays));
output_PM_Nrays = convolution(PM_input,real(CIR_Nrays));
output_TH_Nrays = convolution(TH_input,real(CIR_Nrays));

Ny = size(output_AM_Nrays,1);
ty = linspace(0,Ny,Ny)'*dt;
Npts = 500; 

power_corr_AM_LOS = power_correlation(0.5*abs(output_AM_LOS).^2,Npts);
power_corr_PM_LOS = power_correlation(0.5*abs(output_PM_LOS).^2,Npts);
power_corr_TH_LOS = power_correlation(0.5*abs(output_TH_LOS).^2,Npts);

power_corr_AM_Nrays = power_correlation(0.5*abs(output_AM_Nrays).^2,Npts);
power_corr_PM_Nrays = power_correlation(0.5*abs(output_PM_Nrays).^2,Npts);
power_corr_TH_Nrays = power_correlation(0.5*abs(output_TH_Nrays).^2,Npts);

%% Channel impulse response
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
plot(time,20*log10(abs(CIR_LOS)));
title('Channel Impulse Response for LOS')
xlabel('Time')
ylabel('CIR, dB')
grid on

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
plot(time,20*log10(abs(CIR_Nrays)));
title(strcat('Channel Impulse Response for ',{' '},num2str(Nrays),' rays'))
xlabel('Time')
ylabel('CIR, dB')
grid on

%% Channel transfer function
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
plot(freq/1e9,20*log10(abs(TF_LOS)));
title('Channel Transfer Function')
xlabel('Frequency,GHz')
ylabel('Transfer Function, dB')
grid on

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
plot(freq/1e9,20*log10(abs(TF_Nrays)));
title(strcat('Channel Transfer Function for ',{' '},num2str(Nrays),' rays'))
xlabel('Frequency,GHz')
ylabel('Transfer Function, dB')
grid on

%% Phase-modulated signal
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
subplot(3,1,1)
plot(ts,abs(PM_input))
title('PM Input')
ax = gca;
xlim([0 20e-9])

subplot(3,1,2)
ax = gca;
plot(ty,abs(output_PM_LOS))
title('PM Filter Output (LOS)')
xlabel('Time (sec)')
xlim([0 20e-9])

subplot(3,1,3)
plot(ty,real(output_PM_Nrays))
title('PM Filter Output (N-Rays)')
xlabel('Time (sec)')
xlim([0 20e-9])

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');

subplot(2,1,1)
plot(ts(1:Npts),power_corr_PM_LOS)
title('PM Power Correlation (L.O.S.)')
xlabel('Time Delay (sec)')
ylim([0.95 max(power_corr_PM_LOS)+max(power_corr_PM_LOS)/10])
grid on

subplot(2,1,2)
plot(ts(1:Npts),power_corr_PM_Nrays)
title('PM Power Correlation (N-Rays)')
xlabel('Time Delay (sec)')
ylim([0.95 max(power_corr_PM_Nrays)+max(power_corr_PM_Nrays)/10])
grid on

%% Amplitude-modulated signal
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
subplot(3,1,1)
plot(ts,abs(AM_input))
title('AM Input')
xlim([0 20e-9])

subplot(3,1,2)
plot(ty,abs(output_AM_LOS))
title('AM Filter Output (LOS)')
xlabel('Time (sec)')
xlim([0 20e-9])

subplot(3,1,3)
plot(ty,real(output_AM_Nrays))
title('AM Filter Output (N-Rays)')
xlabel('Time (sec)')
xlim([0 20e-9])

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');

subplot(2,1,1)
plot(ts(1:Npts),power_corr_AM_LOS)
title('AM Power Correlation (L.O.S.)')
xlabel('Time Delay (sec)')
grid on

subplot(2,1,2)
plot(ts(1:Npts),power_corr_AM_Nrays)
title('AM Power Correlation (N-Rays)')
xlabel('Time Delay (sec)')
grid on

%% Thermal signal
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');
subplot(3,1,1)
plot(ts,abs(TH_input))
title('Thermal Input')
xlim([0 20e-9])

subplot(3,1,2)
plot(ty,abs(output_TH_LOS))
title('Thermal Filter Output (LOS)')
xlabel('Time (sec)')
xlim([0 20e-9])

subplot(3,1,3)
plot(ty,real(output_TH_Nrays))
title('Thermal Filter Output (N-Rays)')
xlabel('Time (sec)')
xlim([0 20e-9])

figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', 'on');

subplot(2,1,1)
plot(ts(1:Npts),power_corr_TH_LOS)
title('Thermal Power Correlation (L.O.S.)')
xlabel('Time Delay (sec)')
grid on

subplot(2,1,2)
plot(ts(1:Npts),power_corr_TH_Nrays)
title('Thermal Power Correlation (N-Rays)')
xlabel('Time Delay (sec)')
grid on

function y = IFT(Y,omega,time)
    Nt = size(time,1);
    Nf = size(omega,1);
    y = zeros(Nt,1);
    for i=1:Nt
        y(i) = sum(Y.*exp(1i*omega*time(i))); 
    end
    y = y / Nf; %this is equivalent to having a rectangular window
end

function Y = FT(y,omega,time,sample_rate_multiplier)
    Nf = size(omega,1);
    Y = zeros(Nf,1);
    for i=1:Nf
        Y(i)=sum(y.*exp(-1i*time*omega(i)));
    end
    Y = Y / sample_rate_multiplier;
end

function c = power_correlation(P,Npts)
    c = zeros(Npts,1);
    Pm = mean(P);
    Pmsqr = Pm^2;
    Np = size(P,1);
    for i=1:Npts
        Q = zeros(Np,1);
        Q(i+1:end)=P(1:end-i);
%       Q = circshift(P,i-1);
       Qm = mean(Q);
       c(i) = mean(P.*Q)/Qm/Pm; 
    end
end

function y = convolution(x1,x2)
    conv_num = 3;
    conv_option = ["full","same","valid"];
    y = conv(real(x1), real(x2), conv_option(conv_num)) - ...
        conv(imag(x1), imag(x2), conv_option(conv_num)) + ...
    1i*(conv(real(x1), imag(x2), conv_option(conv_num)) + ...
        conv(imag(x1), real(x2), conv_option(conv_num)));    
end