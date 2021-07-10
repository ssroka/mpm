 %% This script will calculate the radius and temperature evolution of a single saltwater drop

%% ----------------------------- BEGIN ------------------------------------

%%  initialize air and sea meterological quantities

eval(inputFile)
fprintf('\nInitial Conditions Set According to: %s \n',inputFile)
fprintf('%s\n\n',(repmat('-',1,38+length(inputFile))))

if exist('microphysicalConstants.mat','file')
    fprintf('\n Microphysical Constants Set\n')
else
    error('\n Please Set Microphysical Constants\n')
end
fprintf('%s\n\n',(repmat('-',1,38+length(inputFile))))

%% loop through initial radii

% r_fail = false(length(r_0_vec),1);
for r_0_ind = 1:length(r_0_vec)
    clearvars -except inputFile print_flag r_fail r_0_ind Nayar_flag
    eval(inputFile)
    
    r_0 = r_0_vec(r_0_ind);
    fprintf('\n r_0 = %3.2f %sm %s \n',r_0*1e6,char(181),(repmat('-',1,38+length(inputFile))))
    
    % -------------------------------------------------------
    %             Sea Spray Generating Function
    % -------------------------------------------------------
    switch SSGF_str
        case 'singleDrop'
            dFdr0 = 1/r_0*r_0; % 1/ (m^2 s^1 ) units
            Vol_flux = (4.*pi.*r_0.^3/3)*dFdr0;
        case 'Fairall_94'
            % dFdr0 = 3.8*1e-6*U10^(3.4)*5.0*1e-6% Andreas 08
            dFdr0  = SSGF_Fairall94(r_0*1e6);% argument comes in as micrometers
            %         dFdr0_vec(r_0_ind)=dFdr0;
            W_u = 3.8*1e-6*U10^(3.4); % whitecap fraction
            Vol_flux = W_u*dFdr0;
        case 'Ortiz_Suslow_16'
            dFdr0  = SSGF_OrtizSuslow(r_0*1e6);% argument comes in as micrometers
            %         dFdr0_vec(r_0_ind)=dFdr0;
            Vol_flux = 4/3*pi*r_0^3*dFdr0; % units are m^2/(s * micrometers)
    end
    
    % -------------------------------------------------------
    %             Calculate Micro-physical Endpoints
    % -------------------------------------------------------
    helpingAnonFxns;
    s0    = S0/1000; % kg/kg
    T_s   = T_s_0;  % deg Celsius
    [m_s] = compute_initial_drop_conditions(T_s_0,S0,r_0,p0);
    
    storeICs;
    
    try
        compute_Teq;
    catch
	disp('**** Was not able to calculate Teq ****')
        print_comp = true;
    end

    
    try
    compute_req;
    catch
	    disp('**** Was not able to calculate req ****')
    print_comp = true;
    end
    
    try
    compute_tauT;
    catch
	    disp('**** Was not able to calculate tauT ****')
    print_comp = true;
    end
    
    try
    compute_taur;
    catch
	    disp('**** Was not able to calcuate taur ****')
    print_comp = true;
    end
    
    integrate_r_T;
    
    try
    collect_mm_values;
    print_comp = true;
    catch
	    print_comp = false;
	    disp('**** Was not able to collect approximate values ****')
	    T_eq_fmm = NaN;
	    r_eq_fmm = NaN;
	    tau_T_fmm = NaN;
	    tau_r_fmm = NaN;
    end 
    
    tau_f = compute_tauf(U10,T_s_0,r_0,m_s,s0,p0,T_a,maxEr_uf,maxIt);
    
    if print_comp
    print_comparisons;
    end
    saveScript;
    
end
