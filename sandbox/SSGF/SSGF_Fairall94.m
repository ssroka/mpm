%Figure 3 from Failarll et al. 1994, calculating S_v
function [fn_m2] = SSGF_Fairall94(r)
% takes in r in micrometers
% returns the value of Dfdr

% we subdivide Fig 3 into four regions

% the output is in m/(s* micrometer)

% we recall that the points in Fig 3 are in micometers^-2 s^-1 cm^-2
% tracing the  black solid line
FairallWebPlotDig_Fig3;
r0_endpoints = Fairall_Fig3_WebPlotDig_data(:,1);
Sv_endpoints = Fairall_Fig3_WebPlotDig_data(:,2);

Sv = integrate_logspace(r0_endpoints,Sv_endpoints);

try
i_star = find(r0_endpoints<r,1,'last');%index of the lower bound
slope  = (log10(Sv_endpoints(i_star+1))-log10(Sv_endpoints(i_star)))/(log10(r0_endpoints(i_star+1))-log10(r0_endpoints(i_star)));
fn_m2  = 1e-14 * 10.^(log10(Sv_endpoints(i_star))+slope*(log10(r)-log10(r0_endpoints(i_star))));
catch
	keyboard
end
