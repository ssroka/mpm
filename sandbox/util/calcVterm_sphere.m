%% Calculate terminal velocity for a sphere
% [v_term] = calcVterm_sphere(r,rho,g)
% input : 
% r   = radius [m]
% --- optional ---
% rho = density of water [kg/length^3] 
% g   = gravitational acceleration [m/s^2]
%
% output :
% v_term = treminal velocity in air (rho_air = 1 kg/m^3]
%
function [v_term] = calcVterm_sphere(r,rho_w,g)
C = 1; % drag coefficient for Re of 100
rho_a = 1; % kg m^-3 air density
if nargin==1
    rho_w = 1000; % kg m^-3 water density
    g = 9.81;   % m s^-2
elseif nargin==2
    g = 9.81;   % m s^-2
end

m = 4/3*pi*r.^3*rho_w;
A = pi*r.^2;
v_term = sqrt((2.*m*g)./(C*rho_a.*A));

% recalculate Re
nu_air = 2e-5;
Re = v_term*2*r/nu_air;
end
% see
%{
http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450(1969)008%3C0249%3ATVORA%3E2.0.CO%3B2
%}



