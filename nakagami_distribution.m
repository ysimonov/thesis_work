clear all;

z = [0:0.0001:2.5];
Pr = 1;
y = zeros(size(z));
b = zeros(10);
b = ["m = 0.5", "m = 1.0", "m = 1.5", "m = 2.0", ...
    "m = 2.5", "m = 3.0", "m = 3.5", "m = 4.0", ...
    "m = 4.5", "m = 5.0"];

for i = 1:10
m = i / 2;
y = (2 *(m^m) .* z.^(2*m-1)) ./(gamma(m) * Pr^m) .* exp(-m .* z.^2 / Pr);
y(1) = 0;
y(z <= 0) = 0; %additional constraint
plot(z,y,'Marker','.')
hold on
leg = legend(b(1:i))
end 
set(gca,'FontSize',14)
leg.FontSize = 14;
ylabel("f(z)",'FontSize',15)
xlabel("z",'FontSize',15)
grid on
title("Nakagami-m PDF, P_r = 1","FontSize",16)