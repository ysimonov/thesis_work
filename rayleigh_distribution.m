clear all;

z = [0:0.0001:10];
sigma = zeros(size(z));
y = zeros(size(z));
b = zeros(10);
b = ["\sigma = 0.5", "\sigma = 1.0", "\sigma = 1.5", "\sigma = 2.0", ...
    "\sigma = 2.5", "\sigma = 3.0", "\sigma = 3.5", "\sigma = 4.0", ...
    "\sigma = 4.5", "\sigma = 5.0"];

for i = 1:10
sigma = i / 2;
sigma_squarred = sigma.^2;
y = (z ./ sigma_squarred) .* exp(-0.5 * z.^2 / sigma_squarred);
y(z <= 0) = 0; %additional constraint
plot(z,y,'.')
hold on
leg = legend(b(1:i))
end 
set(gca,'FontSize',14)
leg.FontSize = 14;
ylabel("f(z)",'FontSize',15)
xlabel("z",'FontSize',15)
grid on
title("Rayleigh PDF","FontSize",16)