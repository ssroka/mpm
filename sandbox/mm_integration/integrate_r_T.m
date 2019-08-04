% evolve T and r together
%% Another error bound will change anything using the stateVec
% options = odeset('AbsTol',1e-7,'RelTol',2.22045e-7);
%% What is used in plotAndreas05_Fig11
options = odeset('AbsTol',1e-6,'RelTol',2.22045e-6);

P_mb = Pa2mb(p0);
load('microphysicalConstants.mat','R','M_w')
tic
[time_vec, stateVec] = ode45(@(t,stateVec) compute_dsdt(t,stateVec,RH,m_s,T_a,P_mb,t_final,R,M_w,maxEr_s,maxIt),[0 t_final],[r_0;T_s_0],options);
r_t   = stateVec(:,1);
T_s_t = stateVec(:,2);
etime = toc;
fprintf('\n Total Integration Time: %3.1f s\n',etime)
disp(' ')








