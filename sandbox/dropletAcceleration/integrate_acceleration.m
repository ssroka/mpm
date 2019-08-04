% integrate drop acceleration
clear;close all;clc
resultsDir = {'/Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/results/Exp1_01'};
d = dir(sprintf('%s/..',pwd));
for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        addpath(sprintf('%s/../%s',pwd,d(ii).name))
    end
end

n = 0;
for jj = 1:length(resultsDir)
    if ~strcmp(resultsDir{jj}(1:6),'/Users')
        resultsDir{jj} = strcat(pwd,'/',resultsDir{jj});
    end
mydir = dir(resultsDir{jj});
str = {mydir.name};
pat = 'r_0_\w*\.mat';
files_tmp = regexp(str,pat,'match');
for ii = 1:length(files_tmp)
  if ~isempty(files_tmp{ii})
    files{ii+n}{1} = strcat(resultsDir{jj},'/',files_tmp{ii}{1});
  end
end
  n = n + length(files_tmp);
end
for ii = 1:length(files)
    if ~isempty(files{ii})
        load(files{ii}{1})
    
helpingAnonFxns;

t_max = 2; % seconds
u0    = 0; % m/s starting from rest


ra = rho_a(Celsius2Kelvin(ic.T_a),ic.p0);
rs = rho_s(ic.T_s_0,ic.r_0,ic.m_s,ic.S0/1000,ic.p0);
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
[acel_time_vec, u_vec] = ode45(@(t,u) compute_dudt(t,u,ra,rs,ic.U10,ic.r_0,ic.T_a),[0 t_max],u0,options);

tau_ac(ii) = acel_time_vec(max(find(u_vec<ic.U10-exp(-1))));
r_0_vec(ii) = ic.r_0;
plot(acel_time_vec,u_vec)
hold on
    end
end
figure
loglog(r_0_vec*1e6,tau_ac,'o')
set(gcf,'color','w')
set(gca,'xlim',[0.1 10000],'ylim',[1e-6 1e6])
xlabel('r_0 [\mu m]')
ylabel('\tau_{ac} [s]')

for ii = 1:length(d)
    % do not include hidden directories
    if ~strcmp(d(ii).name(1),'.') && d(ii).isdir
        rmpath(sprintf('%s/../%s',pwd,d(ii).name))
    end
end


