
function [dFdr0] = Zhao2006(r0,U10,beta)
% beta is the wave age
if any(r0>500 | r0<30)
    error('Zhao needs 30 < r0 < 500 microns')
end
g = 9.81;
nu = 1.628125802752141e-05;% solved via two curves in fig 7 of paper old value: 1.53e-5;
kappa = 0.4;
u_star_guess = 1;
z0 =@(u_star) 0.0185.*u_star.^2./g;
u_star = fzero(@(u_star) U10-u_star/kappa.*log(10/z0(u_star)),u_star_guess);
CD = 0.5.*U10^(0.5)*1e-3;% (0.8+0.065.*U10)*1e-3; % U10 needs to be > 1 m/s
% beta = g./(w.*U10);
RB = CD.*U10.^3.*beta./(g*nu);
dFdr0 = ...
    7.84e-3*RB^(1.5)*r0.^(-1).*(r0>30).*(r0<=75)+...
    4.41e1*RB^(1.5)*r0.^(-3).*(r0>75).*(r0<=200)+...
    1.41e13*RB^(1.5)*r0.^(-8).*(r0>200).*(r0<=500);


end
