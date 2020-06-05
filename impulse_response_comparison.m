%Author: Yevgeniy Simonov
%Year: 2020

clear all 
close all
experiment = ["dataM/LOS.s2p", ... %EXP(1) Line-of-sight, 3m
              "dataM/NOLOS.s2p", ... %EXP(2) No Line-of-sight, 3m
              "dataM/1R.s2p", ... %EXP(3) 1 Side Reflector, 3m
              "dataM/2R.s2p", ... %EXP(4) 2 Side Reflectors, 3m
              "dataM/3R.s2p", ... %EXP(5) 3 Side Reflectros, 3m
              "dataM/3R_NOLOS.s2p", ... %EXP(6) 3 Side Reflectors No Line-of-sight, 3m
              "dataM/3R_2D.s2p", ... %EXP(7) 3 Side Reflectors, 2 Diffraction Plates, 3m
              "dataM/3R_4D.s2p", ... %EXP(8) 3 Side Reflectors, 4 Diffraction Plates, 3m
              "dataM/4R_4D.s2p"]'; %EXP(9) 4 Side Reflectors, 4 Diffraction Plates, 3m
          
freqs = sparameters("dataM/LOS.s2p").Frequencies;
Nf = size(freqs,1);
Nex = size(experiment,1);

S = {};
for i=1:Nex
    S{end+1} = sparameters(experiment(i)).Parameters;
end

beta = 3.3;
Window = kaiser(Nf,beta);
%Winfreq = freqz(Window,1,Nf,2);
W0 = sum(Window);

TF = zeros(Nf,Nex);
for i=1:Nex
   TF(:,i)= squeeze(S{i}(2,1,:)./(1+S{i}(1,1,:))).*Window; 
end

Nt = Nf;
tstart = 5e-9; %start time
tstop = 50e-9; %end time
time = linspace(tstart, tstop, Nt)';
omega = 2*pi*freqs;

CIR = zeros(Nt,Nex);
for k=1:Nex
    for i=1:Nt
        CIR(i,k) = sum(squeeze(TF(:,k)).*exp(1i*omega*time(i)));
    end
    CIR(:,k) = CIR(:,k) / W0;
end

%% Plot Channel Transfer Function
numexp = (1:9).';
[X,Y] = meshgrid(numexp,freqs);

markers = {'.' '.' '.' '.' '.' '.' '.' '.' '.'};
%markers = {'.' '^' 'diamond' 'p' 'v' 'd' 'x' 'o' '+'};

figure('DefaultAxesFontSize', 14)
gsh = plot3(X,Y/1e9,20*log10(abs(TF))); 
set(gsh,{'Marker'}, markers(:),'MarkerSize',1)
zlabel('20log_{10}|H(f)W(f)|,dB')
ylabel('Frequency, GHz')
xlabel('Experiment Number')
title("Amplitude of Windowed Transfer Function, Kaiser \beta=3.3")
view(140,75)
legend("L.O.S.", ...
       "No L.O.S.", ...
       "1 side reflector", ...
       "2 side reflectors", ...
       "3 side reflectors", ...
       "3 side reflectors, no L.O.S.", ...
       "3 side reflectors, 2 diffractors", ...
       "3 side reflectors, 4 diffractors", ...
       "4 side reflectors, 4 diffractors" );
xlim([1 9])
grid on

%% Plot Channel Impulse Response
[X,Y] = meshgrid(numexp,time);
figure('DefaultAxesFontSize', 14)
gsh = plot3(X,Y*1e9,20*log10(abs(CIR)))
set(gsh,{'Marker'}, markers(:),'MarkerSize',1)
zlabel('20log_{10}|h(t)|,dB')
ylabel('Time, ns')
xlabel('Experiment Number')
title("Amplitudes of Channel Impulse Response, Kaiser \beta=3.3")
view(140,75)
legend("L.O.S.", ...
       "No L.O.S.", ...
       "1 side reflector", ...
       "2 side reflectors", ...
       "3 side reflectors", ...
       "3 side reflectors, no L.O.S.", ...
       "3 side reflectors, 2 diffractors", ...
       "3 side reflectors, 4 diffractors", ...
       "4 side reflectors, 4 diffractors");
zlim([-80 -20])
ylim([5 35])
xlim([1 9])
grid on

%% Identification of multipath waves
CIR_diff = zeros(Nt,Nex-1);
for i=2:Nex
   CIR_diff(:,i-1) = squeeze(sign(abs(CIR(:,i))-abs(CIR(:,1)))).*squeeze(20*log10(smooth(smooth(abs(CIR(:,1))-abs(CIR(:,i)))))); 
end

markers = {'.' '.' '.' '.' '.' '.' '.' '.'};
[X,Y] = meshgrid(numexp(1:end-1),time);
figure('DefaultAxesFontSize', 14)
%set(gca,'xscale','log')
gsh = plot3(X,Y*1e9,CIR_diff)
set(gsh,{'Marker'}, markers(:),'MarkerSize',1)
zlabel('sgn(\Delta h(t))(20log_{10}|h_{n}(t)-h_{1}(t)|),dB')
ylabel('Time, ns')
%xlabel('Experiment Number')
title("Difference Between L.O.S. and Amplitudes of Channel Impulse Response, Kaiser \beta=3.3")
view(140,75)
legend("No L.O.S.", ...
       "1 side reflector", ...
       "2 side reflectors", ...
       "3 side reflectors", ...
       "3 side reflectors, no L.O.S.", ...
       "3 side reflectors, 2 diffractors", ...
       "3 side reflectors, 4 diffractors", ...
       "4 side reflectors, 4 diffractors");
%zlim([-80 -20])
ylim([5 35])
xlim([1 8])
grid on