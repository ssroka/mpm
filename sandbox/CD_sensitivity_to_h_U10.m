
rho_s = 998; % kg/m^3
r0    = 100e-6; % m

mass = 4/3*pi*r0.^3*rho_s;

g = 9.81; % m/s
h = linspace(1,25); % m
U10 = linspace(10,60); % m/s

[GPE,KE] = meshgrid(h*g+U10.^2);
surf(GPE,KE,(GPE+KE))


















