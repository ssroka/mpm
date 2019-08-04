% --- compute r_eq ---------------------------------------------------
%
% reference: Andreas 2005, eq 5.3
%
% input:
%      f = fractional RH (see in.m) RH/100
%      y = models how curvature (the Kelvin effect) and dissolved salt
%      affect the saturation vapor pressure at the surface of an aqueous
%      solution droplet.
%      r_0 = initial droplet radius [m]                 **[in.m]
%      tau_r = time between ejection and the droplet reaching r_eq [m]
%      Andreas 2005, section 5 evaluates sigma_s and rho_s at T_eq
% output:
%  r = the instantaneous droplet radius [m]
%
load('microphysicalConstants.mat','M_w','R','nu','M_s')
f = RH/100; %[none]

% I think there might be a typo in Andreas 2005 eq 5.4 - the second term
% should have a squared denominator

g_req =@(y) (f-1)-y;
dg_reqdr=@(r,sigma_s,rho_s,Phi_s) (2*M_w*sigma_s)/(R*T_a_Kelvin*r^2*rho_w(T_eq-273.15,p0)) - (4*M_w*Phi_s*m_s*nu*r^2*rho_s*pi)/(M_s*(m_s - (4*r^3*rho_s*pi)/3)^2);

if ~exist('r_eq','var')
    r_eq = r_0*2/3; % first iteration
end
m_w_req = compute_m_w(T_eq-273.15,r_eq,m_s,p0,maxEr_s,maxIt);
m_eq = m_s/(M_s*m_w_req);
error_req = abs(g_req(y(T_eq-273.15,T_a,r_eq,m_s,s0,p0)));
r_eq_old = r_eq;
ItNum = 1;
r_eq_vec = zeros(maxIt,1);
r_eq_vec(1) = r_eq_old;

while error_req > maxEr_req && ItNum < maxIt
    m_w_req = compute_m_w(T_eq-273.15,r_eq,m_s,p0,maxEr_s,maxIt);
    s_req = m_s/(m_w_req + m_s);
    m_eq = m_s/(M_s*m_w_req);
    r_eq = r_eq_old-g_req(y(T_eq-273.15,T_a,r_eq,m_s,s_req,p0))/dg_reqdr(r_eq_old,sigma_s(T_eq-273.15,m_s,m_w_req),rho_s(T_eq-273.15,r_eq,m_s,s_req,p0),Phi_s(m_eq));
    error_req = abs(g_req(y(T_eq-273.15,T_a,r_eq,m_s,s_req,p0)));
    r_eq_old = r_eq;
    ItNum = ItNum + 1;
    % look for convergence
    r_eq_vec(ItNum) = r_eq;
end
clear M_w R nu M_s


%% look for convergence
%{
figure
plot(r_eq_vec(r_eq_vec~=0),'*-')
xlabel('iteration number')
ylabel('r_{eq}')
%}










