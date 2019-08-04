function [Q_s, Q_L,time_vec] = integrate_Qs_QL(resultsFiles)



%% ------------------------------------------------------
%            Load in the radius and temperature time histories

load(resultsFiles)

%% ----------------------------- setup paths --------------------------
[pathstr,name,ext] = fileparts(mfilename('fullpath')) ;
added_path_flag = 0;
path_to_fxns = which('rho_s.m');
if isempty(path_to_fxns)
    fprintf('Using files in %s\n',path_to_fxns(1:end-7))
    d = dir(pathstr);
    for ii = 1:length(d)
        % do not include hidden directories
        if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
            addpath(sprintf('%s/%s',pathstr,d(ii).name))
        end
    end
    addpath(pathstr)
    added_path_flag = 1;
end

%% ----------------------------- Begin --------------------------
helpingAnonFxns;

Q_s = zeros(length(time_vec)-1,1);
Q_L = zeros(length(time_vec)-1,1);
m_w = zeros(2,1);
s   = zeros(2,1);
L_v = zeros(2,1);
Volume = zeros(2,1);

for t = 1:length(time_vec)-1
    % Step 1: compute average saltiness [kg/kg]
    m_w(1) = compute_m_w(T_s_t(t),r_t(t),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
    s(1) = ic.m_s/(m_w(1)+ic.m_s);
    m_w(2) = compute_m_w(T_s_t(t+1),r_t(t+1),ic.m_s,ic.p0,ic.maxEr_s,ic.maxIt);
    s(2) = ic.m_s/(m_w(2)+ic.m_s);
    
    % Step 2: compute average volume [m^3]
    Volume(1) = 4/3*pi*r_t(t)^3;
    Volume(2) = 4/3*pi*r_t(t+1)^3;
    Volume_avg = mean(Volume);
    
    % Step 3: compute specific heat [J kg^-1 K^-1]
    cp_avg = 0.5*(SW_SpcHeat(T_s_t(t+1),'C',s(2),'w',ic.p0,'Pa')+SW_SpcHeat(T_s_t(t),'C',s(1),'w',ic.p0,'Pa'));
    
    % Step 4: compute average density over this time [kg m^-3]
    rho_s_avg = 0.5*(rho_s(T_s_t(t+1),r_t(t+1),ic.m_s,s(2),ic.p0) + rho_s(T_s_t(t),r_t(t),ic.m_s,s(1),ic.p0));
    
    % Step 5: compute Q_s
    Q_s(t) = rho_s_avg*Volume_avg*cp_avg*(T_s_t(t)-T_s_t(t+1));
    
    % Step 6: compute average L_v [J/kg]
    L_v(1) = SW_LatentHeat(T_s_t(t),'C',s(1),'w');
    L_v(2) = SW_LatentHeat(T_s_t(t+1),'C',s(2),'w');
    L_v_avg = mean(L_v);
    
    % Step 7: compute Q_L
    Q_L(t) = -rho_s_avg*4/3*pi*r_t(t)^3*(1-(r_t(t+1)/r_t(t))^3)*L_v_avg;
    
    
end

figure
subplot(1,2,1)
yyaxis left
semilogx(time_vec(1:end-1),Q_s);hold on
ylabel('Q_s [J/s]')
yyaxis right
semilogx(time_vec(1:end-1),Q_L);hold on
ylabel('Q_L [J/s]')
xlabel('t [s]')
subplot(1,2,2)
semilogx(time_vec(1:end-1),Q_s+Q_L);hold on

figure(1)
subplot(1,2,1)
yyaxis left
semilogx(time_vec(1:end-1),cumsum(Q_s.*diff(time_vec)./time_vec(2:end)));hold on
ylabel('Integrated Q_s [J]')
yyaxis right
semilogx(time_vec(1:end-1),cumsum(Q_L.*diff(time_vec)./time_vec(2:end)));hold on
ylabel('Integrated Q_L [J]')
xlabel('t [s]')
subplot(1,2,2)
semilogx(time_vec(1:end-1),cumsum(Q_s.*diff(time_vec)./time_vec(2:end))+cumsum(Q_L.*diff(time_vec)./time_vec(2:end)),'displayname',sprintf('r_0=%d %s m',ic.r_0*1e6,char(956)));hold on
legend('-dynamiclegend')

figure(2)
subplot(1,2,1)
yyaxis left
semilogx(time_vec(1:end-1),cumsum(Q_s));hold on
ylabel('Integrated Q_s [J]')
yyaxis right
semilogx(time_vec(1:end-1),cumsum(Q_L));hold on
ylabel('Integrated Q_L [J]')
xlabel('t [s]')
subplot(1,2,2)
semilogx(time_vec(1:end-1),cumsum(Q_s)+cumsum(Q_L),'displayname',sprintf('r_0=%d %s m',ic.r_0*1e6,char(956)));hold on
legend('-dynamiclegend')



figure(3)
yyaxis left
semilogx(time_vec,T_s_t);
ylabel('T')
yyaxis right
semilogx(time_vec,r_t);
ylabel('r')
xlabel('t [s]')



%% ----------------------------- remove user-added paths --------------------------
if added_path_flag
    rm_ASF_paths
end

end






































