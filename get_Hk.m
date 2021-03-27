%clear;close all;clc
function [] = get_Hk(tof_ind)
%{

Author: Sydney Sroka

This script calculates the enthalpy flux and spray CK for a specified set
of environmental conditions
              > relative humidity RH [%]
              > air-sea temperature differnce DT [K]
              > 10-m wind speed U10 [m/s]
Usage:
The user can specify the input parameters in the USER INPUT section below
and the resultant spray-induced enthalpy flux and spra-CK will be
calculated.

%}

%% USER INPUT

% fluid parameters
% rho_w = 1020;    % saltwater density [kg m^-3]
% s0    = 0.034;   % salinity in [kg/kg]
% cp    = 4000;    % specific heat capacity of saltwater [J kg^-1 K^-1]
% Lv    = 2434054; % latent heat of vaporization [J /kg]
% R     = 8.31447; % universal gas constant [J mol^-1 K^-1]
Nayar_flag = true;

% parameters for time-of-flight Newton-Raphson calculation
maxEr_uf = 4e-4; % maximum error in terminal velocity (uf) in [m/s]
maxIt    = 1000; % maximum number of iterations to execute to calc uf

maxEr_s  = 1e-6; % maximum error in salinity (s) in [kg/kg]

% parameters common to all environmental conditions
r0_vec = [50:50:2000];  % initial drop radii [micron]
% SST_C    = 29; % sea surface temperature in [deg C]
% p0       = 100000; % surface pressure [Pa]
% u0 = 0;


% parameters corresponding to SGF from Troitskaya et al. 2018a
eq_N = 9; % eq_N \in [5,7,9] - which equation to use for the number of bags
% 5 = derived equation
% 7 = laboratory conditions
% 9 = field conditions
Omega = 2.5; % wave-age parameter (usually between 2.5- 3.5)

% environmental conditions

% _vec quantities will be looped over holding the other parameters
% constant at their _default values
U10_vec = [20:10:60];    % [m/s]
RH_vec  = [88:2:98];    % [%]
DT_vec  = [0.5:0.5:3.5];  % [K]
SST_vec  = [27 28 29];  % [K]


% approximation flag
% false => use the full microphysical model formula
%          WARNING if set to false and using results in SST_C_29
%                  r0 MUST be \in [50:50:1200]
%                  DT MUST be \in [1:1:3]
%                  RH MUST be \in [80:5:95]
%          WARNING if set to false and using results in SST_C_27
%                  r0 MUST be \in [100:50:1200]
%                  DT MUST be \in [0.5:0.5:3.5]
%                  RH MUST be \in [88:2:96]
%
% true => use the fitted equation
use_approx_flag = false;

r0_vec_sort = [31 2 13 24 25 26 27 28 29 30 32 33 34 35 36 37 38 39 40 1 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23];
r0_vec_ALL = [50:50:2000];

tof_vec = [0.1:0.1:2];
%% ---------------- Begin -------------------------------------------------

action_str = 'add';
paths_for_calc_CK;

helpingAnonFxns;

nDT = length(DT_vec);
nr0 = length(r0_vec);
nRH = length(RH_vec);
nU10 = length(U10_vec);


% parameters that are fixed given the directory
% arbitrarily select the first values in r0_vec, DT_vec, and RH_vec
TH_src = sprintf('/home/ssroka/mpm/sandbox/results/cp_of_timehistory_201_thru_466/SST_27/timehistory_%dK_%d_RH',DT_vec(1)*10,RH_vec(1));
load(TH_src,'timehistory'); % will instantiate a variable called timehistory

s0 = timehistory(1).ic.s0;        % [kg/kg]
p0 = timehistory(1).ic.p0;        % [Pa]


ZERO = zeros(nDT,nRH,nU10);


    tof_fac = tof_vec(tof_ind);
    
    for SST_ind = 1:length(SST_vec)
        SST_C = SST_vec(SST_ind);
        % microphysical model data is stored
        TH_src = sprintf('/home/ssroka/mpm/sandbox/results/cp_of_timehistory_201_thru_466/SST_%d',SST_C);
        
        HK = ZERO;
        HK_Troit = ZERO;
        
        for DT_ind = 1:nDT
            DT = DT_vec(DT_ind);
            fprintf('DT = %4.1f  -----------------------\n',DT)
            
            % parameters that are fixed given DT
            % arbitrarily select the first values in r0_vec and RH_vec
            TH = sprintf('%s/timehistory_%dK_%d_RH',TH_src,DT*10,RH_vec(1));
            load(TH,'timehistory'); % will instantiate a variable called timehistory
            
            T_a = timehistory(1).ic.T_a;      % [deg C]
            
            for RH_ind = 1:nRH
                RH = RH_vec(RH_ind);
                fprintf('RH = %4.1f\n\n',RH)

                % parameters that are fixed given RH
                TH = sprintf('%s/timehistory_%dK_%d_RH',TH_src,DT*10,RH);
                load(TH,'timehistory'); % will instantiate a variable called timehistory
                
                for U10_ind = nU10:-1:1
                    U10 = U10_vec(U10_ind);

                    dFdr_vec = calc_Troit_SGF(U10,eq_N,Omega,r0_vec*1e-6); % meters
                                        
                    for r0_ind = 1:nr0
                        
                        r0_file = r0_vec_sort(round(r0_vec(r0_ind))==r0_vec_ALL);
                        
                        r0_m = r0_vec(r0_ind)*1e-6;              % convert to [m]
                        
                        m_s = timehistory(r0_file).ic.m_s;      % [kg]
                        t_full = timehistory(r0_file).time_vec;
                        r_full = timehistory(r0_file).r_t;
                        T_full = timehistory(r0_file).T_s_t;
                        
                        tof = tof_fac.*compute_tauf(U10,SST_C,r0_m,m_s,s0,p0,T_a,maxEr_uf,maxIt); % meters
                        inds = t_full<tof;
                        
                        Tfinal = interp1(t_full,T_full,tof);
                        rfinal = interp1(t_full,r_full,tof);
                        t = [t_full(inds); tof];
                        T = [T_full(inds); Tfinal];
                        r = [r_full(inds); rfinal];
                        
                        T_K = Celsius2Kelvin(T);
                        
                        clear t_full r_full T_full
                        
                        
                        rho_0 = SW_Density(T(1),'C',s0,'w',p0,'Pa'); % initial saltwater drop density
                        m0 = 4/3*pi*r0_m.^3*rho_0;
                        
                        rho_w = SW_Density(0.5*(T(1)+T(end)),'C',0,'w',p0,'Pa');  % average freshwater drop density
                        m_v = 4/3*pi*(r0_m.^3-r(end).^3)*rho_w;   % mass of the drop that became vapor
                        
                        m_f = m0 - m_v;  % final mass
                        
                        cpw = SW_SpcHeat(T(1),'C',s0,'w',p0,'Pa'); % the specific heat of the initial satlwater drop
                        
                        % this is multiplying water vapor mass which
                        % doesn't have any salt content
                        Lv0 = SW_LatentHeat(T(1),'C',0,'w');
                        Lvf = SW_LatentHeat(T(end),'C',0,'w');
                        
                        %                         HK0 = m0*T_K(1)*cpw; % initial Hk
                        
                        %                         HKf = m_f*T_K(end)*cpw + m_v*Lvf + m_v*T_K(end)*cpw; % final Hk
                        
                        
                        HK_drop(r0_ind) = (m0*T_K(1)*cpw - m_f*T_K(end)*cpw) + m_v*Lvf + m_v*T_K(end)*cpw;
                        
                        m_troit = SW_Density(T_K(1),'C',s0,'w',p0,'Pa').*4./3.*pi.*r.^3;
                        QK_Troit(r0_ind) = 4018*(m_troit(1)*(T_K(1)-min(T_K))+m_troit(end)*(min(T_K)-T_K(end))); % 32

                    end
                    
                    HK(DT_ind,RH_ind,U10_ind) = trapz(r0_vec,HK_drop.*dFdr_vec);
                    
                    HK_Troit(DT_ind,RH_ind,U10_ind) = trapz(r0_vec,QK_Troit.*dFdr_vec);
        
                end
            end
        end
        save(sprintf('HK_SST_%d_tf_%d',SST_C,tof_ind))
        clc
    end



action_str = 'remove';
paths_for_calc_CK;

% figure(1)
% update_figure_paper_size
% print('HK_Troit','-dpdf')
