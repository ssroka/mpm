
function [m_w, s] = compute_m_w(T,r,m_s,P,maxEr_s,maxIt)
% function [m_w] = compute_m_w(T,r,m_s)
% compute m_w in kg from r [m] and T [degrees C]
%
% A system of 3 equations and 3 unknowns
% M_total = m_s + m_w
% M_total = 4 pi r^3 / 3 * rho_s
% rho_s = f(m_w,m_s,M_s,rho_w,v_a(T,r))
% results in a quadratic equation for m_w
% m_w^2 + [m_s + v_a(T,r) rho_w m_s /M_s - rho_w 4 pi r^3 /3 ] m_w - rho_w 4 pi r^3 m_s /3
%
% input:
%     T   = water temperature    [deg C]
%     r   = radius of drop [m]
%     m_s = mass of salt in drop [kg]
%     P   = pressure in [Pa]
%
% constants:
%    M_s  = molar mass of salt [kg/mol]
%
% output:
%    m_w = mass of pure water  [kg]
%

global Nayar_flag

load('microphysicalConstants.mat','M_s')

if Nayar_flag
    % Use Andreas 1989 functions to compute starting salinity
    quadFormula =@(a,b,c) (-[b;b] + [sqrt(b.^2-4.*a.*c);-sqrt(b.^2-4.*a.*c)])./(2.*[a;a]);
    
    % rho_w
%     rho_w = (999.8396 + 18.224944 * T- 7.922210 * 10^(-3)* T.^2)./(1+1.8159725*10^(-2).*T); %eq 13a
    rho_w = SW_Density(T,'C',0*T,'w',P,'Pa');

    % v_a(T,r)
    v_a0 =@(T)   1e-6*(12.97+0.2340*T-4.210*10^(-3)*T.^2 + 2.857*10^-5 * T.^3); %eq 16
    Sv   =@(T)   10^-6*( 2.982 - 4.970 * 10^(-2)*T + 6.032 * 10^(-4)*T.^2); %eq 17
    c    =@(r)   10^(-3) * m_s./M_s./(4*pi*r.^3/3); % eq 15
    v_a  =@(T,r) v_a0(T) + Sv(T).*sqrt(c(r)); % eq 14
    
    m_w_a = ones(length(r),1);
    m_w_b = m_s + v_a(T,r) .* rho_w .* m_s ./M_s - rho_w .* 4 .* pi .* r.^3 ./3 ;
    m_w_c = - rho_w .* 4 .* pi .* r.^3 .* m_s ./3 ;
    
    m_w_vec = quadFormula(m_w_a,m_w_b,m_w_c);
    if sum(double(m_w_vec>0))==2
        error('both m_w values are positive')
    end
    m_w = m_w_vec(m_w_vec>0);
    s_old = m_s/(m_w+m_s);
    
    
    % use Newton's Method to solve for new salinity
    Iter = 1;
    if isempty(s_old)
	error('compute m_w, s_old is empty')
    end
    rho_s_old_tmp = SW_Density(T,'C',s_old,'w',P,'Pa');
    err = m_s/s_old-4/3*pi*r^3*rho_s_old_tmp;
        
    
    
    while Iter< maxIt && err > maxEr_s
        if length(s_old)~=length(T) || length(P)~=length(T)
	error('length mismatch between s and T')
        end
        rho_s_old = SW_Density(T,'C',s_old,'w',P,'Pa');
        f = m_s/s_old-4/3*pi*r^3*rho_s_old;
        
        d_rho_d_s = SW_Density_ds(T,'C',s_old,'w',P,'Pa');
        dfds = -m_s/(s_old.^2)-4/3*pi*r^3*d_rho_d_s;
        
        s_new = s_old - f/dfds;
        
        rho_s_tmp = SW_Density_ds(T,'C',s_new,'w',P,'Pa');
        f_tmp = m_s/s_new-4/3*pi*r^3*rho_s_tmp;
        
        err = abs(f_tmp);
        
        s_old = s_new;
        Iter = Iter + 1;
    end
    
    m_w = m_s*(1 - s_old)/(s_old);
    s = s_old;
else
    
    quadFormula =@(a,b,c) (-[b;b] + [sqrt(b.^2-4.*a.*c);-sqrt(b.^2-4.*a.*c)])./(2.*[a;a]);
    
    % rho_w
    rho_w = (999.8396 + 18.224944 * T- 7.922210 * 10^(-3)* T.^2)./(1+1.8159725*10^(-2).*T); %eq 13a
    
    % v_a(T,r)
    v_a0 =@(T)   1e-6*(12.97+0.2340*T-4.210*10^(-3)*T.^2 + 2.857*10^-5 * T.^3); %eq 16
    Sv   =@(T)   10^-6*( 2.982 - 4.970 * 10^(-2)*T + 6.032 * 10^(-4)*T.^2); %eq 17
    c    =@(r)   10^(-3) * m_s./M_s./(4*pi*r.^3/3); % eq 15
    v_a  =@(T,r) v_a0(T) + Sv(T).*sqrt(c(r)); % eq 14
    
    m_w_a = ones(length(r),1);
    m_w_b = m_s + v_a(T,r) .* rho_w .* m_s ./M_s - rho_w .* 4 .* pi .* r.^3 ./3 ;
    m_w_c = - rho_w .* 4 .* pi .* r.^3 .* m_s ./3 ;
    
    m_w_vec = quadFormula(m_w_a,m_w_b,m_w_c);
    if sum(double(m_w_vec>0))==2
        error('both m_w values are positive')
    end
    m_w = m_w_vec(m_w_vec>0);
    s = m_s/(m_w+m_s);
end
end


