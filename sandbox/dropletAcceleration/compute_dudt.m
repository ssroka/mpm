function [dudt] = compute_dudt(t,u,ra,rs,U10,r_0,T_a)


% ra = density of air, abbreviated from rho_a to avoid recomputation at
%      each time step
% ra = density of water, abbreviated from rho_a to avoid recomputation at
%      each time step


if isnan(u)
	error('u has NaNs')
end
dudt = 3/8 * ra/rs*(U10 - u)^2/r_0*c_drag_clift(U10-u,2*r_0,nu_a(T_a));















end
