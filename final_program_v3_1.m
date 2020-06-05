close all
clear all

savefiles = true;

Nt = 3000; %Nf; %number of time domain points
Npts = 500; %number of delay points for power correlation
time_multiplier = 20; %Duration of signal in relation to impulse response

Nexp = 4; %Select experiment number
experiment = ["dataM/LOS.s2p", ... %EXP(1) Line-of-sight,3m
              "dataM/NOLOS.s2p", ... %EXP(2) No Line-of-sight, 3m
              "dataM/1R.s2p", ... %EXP(3) 1 Side Reflector, 3m
              "dataM/2R.s2p", ... %EXP(4) 2 Side Reflectors, 3m
              "dataM/3R.s2p", ... %EXP(5) 3 Side Reflectros, 3m
              "dataM/3R_NOLOS.s2p", ... %EXP(6) 3 Side Reflectors No Line-of-sight, 3m
              "dataM/3R_2D.s2p", ... %EXP(7) 3 Side Reflectors, 2 Diffraction Plates, 3m
              "dataM/3R_4D.s2p", ... %EXP(8) 3 Side Reflectors, 4 Diffraction Plates, 3m
              "dataM/4R_4D.s2p"]; %EXP(9) 4 Side Reflectors, 4 Diffraction Plates

%display figures if not saving them
if(savefiles)
    visible_opt = 'off';
else
    visible_opt = 'on';
end

power_level = 0; %in dBm (Used in Measurements)

%Cutoff TF with high frequency response deviation
%Can be used to remove phase deviation at the start of frequency response
remove_frequency_deviation = true;
 
S = sparameters(experiment(Nexp)).Parameters;
freq1 = sparameters(experiment(Nexp)).Frequencies;

%scattering two-port parameters
S21 = squeeze(S(2,1,:));

%characteristic impedance
Z0 = 50; 

%channel transfer function
TF1 = S21/2; 

%remove lower frequency components that cause phase deviation
%the upper frequencies have to be removed too to keep  
%center frequency at the same position
freq_dev = 0.5e9; 
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
freq_GHz = freq/1e9;
omega = 2*pi*freq;
Nf = size(freq,1);
    
%the number of times the sampling rate is greater than CIR's bandwidth
df = freq(2)-freq(1); %frequency step
Tspan = 1/df; %maximum time span of the Fourier transform 

%the time domain boundaries were identified by inspection
tmin_arr = [11.06; 6.6; 5.96; 11.24; 11.38; 0; 2.051; 5.203; 0]*1e-9;
tmax_arr = [17.6; 33.38; 26.03; 25.97; 26.05; 31.77; 65.33; 64.58; 53.93]*1e-9;
tmin = tmin_arr(Nexp);
tmax = tmax_arr(Nexp);

time = linspace(tmin,tmax,Nt).';
dt = time(2)-time(1);
time_ns = time*1e9;

%used to correct amplitude of the continuous sampled Fourier transform
sample_rate_multiplier = Tspan / (tmax-tmin) * Nt / Nf; 

%find channel impulse response
CIR = IFT(TF,omega,time);

CIR_dB = 20*log10(abs(CIR));

clear S A B S12 S22 S11 freq1 S21 
%% Plot CIR and Transfer Function

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
plot(time_ns,20*log10(abs(CIR)));
grid on
ylabel('Channel Impulse Response, dB')
xlabel('Time, ns')

save_to_png(gcf,'Impulse_Response',Nexp,savefiles)

markers = {'*' '.'};
gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
gsh = plot(time_ns,real(CIR),time_ns,imag(CIR));
grid on
ylabel('Channel Impulse Response, linear')
xlabel('Time, ns')
set(gsh,{'Marker'}, markers(:),'MarkerSize',4)
legend('Re(CIR)','Im(CIR)')

save_to_png(gcf,'Impulse_Response_Re_Im',Nexp,savefiles)

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
TF_Test = FT(CIR,omega,time,sample_rate_multiplier);
plot(freq_GHz,20*log10(abs(TF_Test)),'.','MarkerSize',20)
hold on
plot(freq_GHz,20*log10(abs(TF)),'x','MarkerSize',2)
grid on
legend('Fourier Transform of CIR','Transfer Function')
xlabel('Frequency, Hz')
ylabel('Transfer Function, dB')

save_to_png(gcf,'Transfer_Function',Nexp,savefiles)

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
subplot(2,1,1)
plot(freq_GHz,angle(TF))
grid on
xlabel('Frequency, GHz')
ylabel('Phase (Deg)')
title('Phase of the Transfer Function')
subplot(2,1,2)
plot(freq_GHz,unwrap(angle(TF)),'x','MarkerSize',2)
grid on
xlabel('Frequency, GHz')
ylabel('Unwrapper Phase (Deg)')

save_to_png(gcf,'Transfer_Function_phase',Nexp,savefiles)

%% GENERATE DIFFERENT SIGNALS AND ANALYZE

%assume that the noise floor is defined by the first data point of CIR
SNR = abs(CIR_dB(1))-power_level; 

fc = (freq(end)-freq(1))/2+freq(1);
wc = 2*pi*fc;

Ns = ceil(Nt*time_multiplier);
ts = linspace(0,Ns,Ns).'*dt;
lowpass_to_bandpass = exp(1i*wc*ts);

%% THERMAL SOURCE
x_th = randn(Ns,1)+cos(wc*ts);
y_th = convolution(x_th,real(CIR));

Ny = size(y_th,1);
ty = linspace(0,Ny,Ny).'*dt;

P_th = 0.5*abs(y_th).^2;
P_in_th = 0.5*abs(x_th).^2;

P_th_corr = power_correlation(P_th, Npts);
P_in_th_corr = power_correlation(P_in_th, Npts);

%% AMPLITUDE MODULATED SIGNAL
x_am = (1+unifrnd(-1,1,Ns,1)).*cos(2*pi*fc*ts);
y_am = convolution(x_am,real(CIR));

P_am = 0.5*abs(y_am).^2;
P_in_am = 0.5*abs(x_am).^2;

P_am_corr = power_correlation(P_am, Npts);
P_in_am_corr = power_correlation(P_in_am, Npts);

%% PHASE MODULATED SIGNAL

phase = unwrap(wc*ts) + unwrap(2*pi*rand(Ns,1));
x_pm = real(exp(1i*phase));
y_pm = convolution(x_pm,real(CIR));

P_pm = 0.5*abs(y_pm).^2;
P_in_pm = 0.5*abs(x_pm).^2;

P_pm_corr = power_correlation(P_pm, Npts);
P_in_pm_corr = power_correlation(P_in_pm, Npts);

%% POWER SPECTRAL DENSITY (in Vrms^2/Hz)
df = freq(2)-freq(1); %required scaling

%Refer to Wiener-Khinchin theorem for more details 
ACF_th = autocorr(y_th/sqrt(2),'NumLags',Ny-1)*var(y_th/sqrt(2),1);
ACF_am = autocorr(y_am/sqrt(2),'NumLags',Ny-1)*var(y_am/sqrt(2),1);
ACF_pm = autocorr(y_pm/sqrt(2),'NumLags',Ny-1)*var(y_pm/sqrt(2),1);

Nacf = size(ACF_th,1);
tau1 = linspace(0,Nacf,Nacf).'*dt;

sample_rate_multiplier1 = Tspan / (tau1(end)-tau1(1)) * Nacf / Nf;

PSD_th = FT(ACF_th,omega,tau1,sample_rate_multiplier1)/df;
PSD_am = FT(ACF_am,omega,tau1,sample_rate_multiplier1)/df;
PSD_pm = FT(ACF_pm,omega,tau1,sample_rate_multiplier1)/df;

%% PLOT SIGNALS (BANDPASS AND LOWPASS EQUIVALENT ENVELOPE)
Nspl = floor(size(ts,1)/10);
tsp = ts(1:Nspl)*1e9;

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);

subplot(3,1,1)
plot(tsp, x_th(1:Nspl),'--')
grid on
xlabel('Time, ns')
ylabel('Thermal Signal, V')
xlim([0 tsp(end)])
title('Transmitted Signals')

subplot(3,1,2)
plot(tsp, x_am(1:Nspl),'--')
grid on
xlabel('Time, ns')
ylabel('AM Signal, V')
xlim([0 tsp(end)])

subplot(3,1,3)
plot(tsp, x_pm(1:Nspl),'--')
grid on
xlabel('Time, ns')
ylabel('PM Signal, V')
xlim([0 tsp(end)])

save_to_png(gcf,'Transmitted_signals',Nexp,savefiles)

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);

subplot(3,1,1)
plot(tsp, y_th(1:Nspl),'--')
grid on
xlabel('Time, ns')
ylabel('Thermal Signal, V')
xlim([0 tsp(end)])
title('Bandpass Representation of Received Signals')

subplot(3,1,2)
plot(tsp, y_am(1:Nspl),'--')
grid on
xlabel('Time, ns')
ylabel('AM Signal, V')
xlim([0 tsp(end)])

subplot(3,1,3)
plot(tsp, y_pm(1:Nspl),'--')
grid on
xlabel('Time, ns')
ylabel('PM Signal, V')
xlim([0 tsp(end)])

save_to_png(gcf,'Received_signals',Nexp,savefiles)

%% PLOTTING ANALYZED DATA
%autocorrelation functions 
gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
plot(tau1*1e9,ACF_am,'.')
hold on
plot(tau1*1e9,ACF_th,'.')
hold on
plot(tau1*1e9,ACF_pm,'*','MarkerSize',4,'color','k')
grid on
legend('Amplitude modulated Signal', ...
       'Thermal Noise', ...
       'Phase modulated Signal')
xlabel('Delay Time, ns')
ylabel('ACF(y)')
hold off

save_to_png(gcf,'Autocorrelation_Voltage',Nexp,savefiles)

%%
%correlation coefficients 
gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);

N1 = size(P_pm_corr,1);
tau = linspace(0,N1,N1).'*dt*1e9; %in ns

plot(tau, P_pm_corr, '-*', 'color', 'm');
hold on
plot(tau, P_am_corr, '-^', 'color', 'red');
hold on
plot(tau, P_th_corr, '-x', 'color', 'b');
grid on
legend('Phase modulated Signal', ...
       'Amplitude modulated Signal', ...
       'Thermal Noise')
xlabel('Delay Time, ns')
ylabel('Power Correlation')
title('Power Correlation at the Receiver')

save_to_png(gcf,'Power_Correlation_Receiver',Nexp,savefiles)

N2 = size(P_in_pm_corr,1);
tau = linspace(0,N2,N2).'*dt*1e9; %in ns

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
plot(tau, P_in_pm_corr, '-x', 'color', 'm');
hold on
plot(tau, P_in_am_corr, '--', 'color', 'r');
hold on
plot(tau, P_in_th_corr, '-.', 'color', 'b');
grid on
legend('Phase modulated Signal', ...
       'Amplitude modulated Signal', ...
       'Thermal Noise')
xlabel('Delay Time, ns')
ylabel('Power Correlation')
title('Power Correlation at the Transmitter')

save_to_png(gcf,'Power_Correlation_Transmitter',Nexp,savefiles)

%%
%power spectral density
gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
plot(freq_GHz, 10*log10(abs(PSD_pm)), '>', ...
     freq_GHz, 10*log10(abs(PSD_am)), '+', ...
     freq_GHz, 10*log10(abs(PSD_th)), 'x', 'MarkerSize', 2)
grid on
legend('Phase modulated Signal', ...
       'Amplitude modulated Signal', ...
       'Thermal Noise')
xlabel('Frequency, GHz')
ylabel(strcat('PSD of the Received Signal (SNR = ',num2str(SNR),'dB) dB(V_{rms}^{2}/Hz)'))

save_to_png(gcf,'Power_Spectral_Density',Nexp,savefiles)

%input power / output power

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
plot(ty*1e9, 10*log10(P_am), '.', ...
     ty*1e9, 10*log10(P_th), '.', ...
     ty*1e9, 10*log10(P_pm), '.') %ts*1e9, P_in_pm, '*', 
xlabel('Time, ns')
ylabel('Normalized Instantaneous Signal Power at the Receiver, dB(W.\Omega)')
legend('AM','Thermal','PM')
title('Power of Different Signals at the Receiver')
grid on
save_to_png(gcf,'Instantaneous_Power',Nexp,savefiles)

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
p = histogram(abs(y_th),'Normalization','probability','BinWidth',0.0008); 
title('Probability Density Function (Thermal Source)')
grid on
xlabel('Envelope Magnitude')

save_to_png(gcf,'PDF_TH',Nexp,savefiles)

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
p = histogram(abs(y_am),'Normalization','probability','BinWidth',0.0008); 
title('Probability Density Function (Amplitude Modulation)')
grid on
xlabel('Envelope Magnitude')

save_to_png(gcf,'PDF_AM',Nexp,savefiles)

gcf = figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1], 'Visible', visible_opt);
p = histogram(abs(y_pm),'Normalization','probability','BinWidth',0.0008); 
title('Probability Density Function (Phase Modulation)')
grid on
xlabel('Envelope Magnitude, V')

save_to_png(gcf,'PDF_PM',Nexp,savefiles)

function Y = FT(y,omega,time,sample_rate_multiplier)
    Nf = size(omega,1);
    Y = zeros(Nf,1);
    for i=1:Nf
        Y(i)=sum(y.*exp(-1i*time*omega(i)));
    end
    Y = Y / sample_rate_multiplier;
end

function y = IFT(Y,omega,time) %for filters
    Nt = size(time,1);
    Nf = size(omega,1);
    y = zeros(Nt,1);
    for i=1:Nt
        y(i) = sum(Y.*exp(1i*omega*time(i))); 
    end
    y = y / Nf; %this is equivalent to having a rectangular window
end

function [y,freqnew,Nfnew] = LOWPASS_IFT(Y,freq,time) %for real networks
    Nt = size(time,1);
    Nf = size(freq,1);
    f1 = freq(1);
    df = freq(2)-f1;
    Tmax = 1/df;
    Nzp = f1/df;
    Ynew = vertcat(zeros(Nzp,1),Y);
    Nfnew = Nf + Nzp;
    freqnew = linspace(0,Nfnew,Nfnew).'*df;
    y = zeros(Nt,1);
    for i=1:Nt
        y(i) = 2*real(sum(Ynew(2:end).*exp(1i*freqnew(2:end)*time(i)))); 
    end
    y = y + Ynew(1);
    y = y / Nfnew; %this is equivalent to having a rectangular window
end

function c = power_correlation(P,Npts)
    c = zeros(Npts,1);
    Pm = mean(P);
    Pmsqr = Pm^2;
    Np = size(P,1);
    for i=1:Npts
       Q = zeros(Np,1);
       Q(i+1:end)=P(1:end-i);
       %Q = circshift(P,i-1);
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

function save_to_png(gcf, name, Nexp, savefiles)
    experiments = {'LOS', ... 
                  'NOLOS', ... 
                  '1REFL_LOS', ... 
                  '2REFL_LOS', ... 
                  '3REFL_LOS', ... 
                  '3REFL_NOLOS', ... 
                  '3REFL_2DIFF_LOS', ... 
                  '3REFL_4DIFF_LOS', ... 
                  '4REFL_4DIFF_LOS'}';
    subdirectory = 'results';
    if ~exist(subdirectory, 'dir')
        [status, msg, msgID] = mkdir(subdirectory)
    end
    subdirectory = strcat('results/',experiments{Nexp});
    if ~exist(subdirectory,'dir')
        [status, msg, msgID] = mkdir(subdirectory)
    end
    if(savefiles)
        path = strcat('/',subdirectory,'/',name);
        saveas(gcf,[pwd path]); %save as matlab figure
        saveas(gcf,[pwd path],'png'); %save as png
    end
end