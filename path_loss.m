clear all
close all

S = {};
%data for path-loss calculations (1001 pts)
freqs_1001pts = sparameters("dataM/100cm_SEP.s2p").Frequencies;
N2 = size(freqs_1001pts,1);
S{end+1} = sparameters("dataM/100cm_SEP.s2p").Parameters;
S{end+1} = sparameters("dataM/194cm_SEP.s2p").Parameters;
S{end+1} = sparameters("dataM/363cm_SEP.s2p").Parameters;

freq_GHz = freqs_1001pts / 1e9;

%calculate LOS gain for path-loss calculations
Gain = zeros(N2,3);
Gain(:,1) = abs(squeeze(S{end-2}(2,1,:)/2)).^2; %1.00m separation
Gain(:,2) = abs(squeeze(S{end-1}(2,1,:)/2)).^2; %1.94m separation
Gain(:,3) = abs(squeeze(S{end}(2,1,:)/2)).^2; %3.63m separation

%calculate path-loss coefficient 
Pt = 1e-3; %transmitted power
Pr = Gain * Pt; %received power

alpha = zeros(N2,2);
d0 = 1.00; %reference distance
d1 = 1.94;
d2 = 3.63; 
alpha(:,1) = smoothdata(smoothdata(log(Pr(:,2)./Pr(:,1))./log(d0/d1)));
alpha(:,2) = smoothdata(smoothdata(log(Pr(:,3)./Pr(:,1))./log(d0/d2)));
alpha_mean = squeeze(alpha(:,1) + alpha(:,2))/2;

d = linspace(1,100)'; %distance in meters
Nd = size(d,1);
Pd = zeros(N2,Nd);

for i=1:Nd
    Pd(:,i) = Pr(:,1).*(d0./d(i)).^alpha_mean(:);
end

[X,Y]=meshgrid(d,freqs_1001pts);
%%
figure('DefaultAxesFontSize', 14, 'units',...
    'normalized','outerposition',[0 0 1 1]);
plot(freq_GHz,alpha(:,1),freq_GHz,alpha(:,2),freq_GHz,alpha_mean)
legend('1.94m','3.63m','Mean')
title('Path-loss Exponent w.r.t. Reference 1m Separation')
ylabel('Path-loss Coefficient')
xlabel('Frequency, GHz')
grid on

figure('DefaultAxesFontSize', 14)
mesh(X,Y/1e9,10*log10(Pd))
colorbar
colormap jet
zlabel('10log_{10}[P(d,f)],dB')
ylabel('Frequency, GHz')
xlabel('Distance, m')
title("Expected Recever Power as a Function of Distance")
view(135,55)

%%