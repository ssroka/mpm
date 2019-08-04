


% Define ranges for each quantity

%% EDIT THESE WHEN RUNNING AGAIN
% TODO: fix pressure

count=1;
% explore r0 variations
% ----------- other parameters --------------------------------------------
% input file for Approximation forumlas for the microphysical properties of
% saline droplets
RH = 90; % relative humidity

% ambient air temperature
T_a = 18; % Celsius

% initial droplet temperature (and local SST)
T_s = 20; % Celsius

% salinity of drop
S = 34; % ppt

% pressure of air
p0 = 100000; % Pa

% initial droplet radius
r_0_vec = 100*1e-6; % m

% 10-m wind speed
U10 = 20; % m/s
% time
t = 0.1; %s

% Newton-Raphson
% for both u_f and r_eq
maxIt = 100; % the maxium number of iterations before aborting the iterative scheme

maxEr_req = max(r_0_vec)*0.0001; % the max error in r_eq
maxEr_uf  = calcVterm_sphere(max(r_0_vec))*0.001;  % the max error in u_f

% micro-physical endpoints from caption in Fig 11 of Andreas 2005
T_eq_exact = 17.07;% deg C
r_eq_exact = 61.44 * 10^-6;% m
tau_T_exact = 0.176;% s
tau_r_exact = 303;% s

%%

for trial = 1:15
    switch trial
        case 1
            Andreas05_Fig11_in
            
            try
                % ------------- vary r_0 ------------------------------------------------------
                for r_0_vec = [10:5:400]*1e-6
                    r_0_vec
                    clearvars -except count r_0_vec Ck_mat RH T_a T_s S p0 U10 Ck_vec1 maxIt maxEr_req maxEr_uf tau_r_exact

                    Ck_vec1(count)=C_k(1);
                    count = count+1;
                end
                save('Ck_mat1','Ck_vec1');
                
            catch
                continue
            end
        case 2
            Andreas05_Fig11_in
            
            try
                % ------------- vary r_0 with S------------------------------------------------------
                for r_0_vec = [10:5:400]*1e-6
                    r_0_vec
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec2 maxIt maxEr_req maxEr_uf tau_r_exact
                    countS = 1;
                    for S = 30:40
                        calcCk_1drop
                        Ck_vec2(count,countS)=C_k(1);
                        countS=countS+1;
                    end
                    % r_0 with S
                    count = count+1;
                end
                save('Ck_mat2','Ck_vec2');
            catch
                continue
            end
        case 3
            Andreas05_Fig11_in
            
            try
                % ------------- vary r_0 with U10 ------------------------------------------------------
                count =1;
                S = 34;%psu
                for r_0_vec = [10:5:400]*1e-6
                    r_0_vec
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec3 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for U10 =10:5:80
                        calcCk_1drop;
                        Ck_vec3(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat3','Ck_vec3');
            catch
                continue
            end
        case 4
            Andreas05_Fig11_in
            
            try
                % ------------- vary r_0 with RH  ------------------------------------------------------
                count =1;
                U10 = 20;% m/s
                for r_0_vec = [10:5:400]*1e-6
                    r_0_vec
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec4 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for RH =70:90
                        calcCk_1drop;
                        Ck_vec4(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat4','Ck_vec4');
            catch
                continue
            end
        case 5
            Andreas05_Fig11_in
            
            try
                % ------------- vary r_0 with T_s  ------------------------------------------------------
                count =1;
                RH = 90; %percent
                for r_0_vec = [10:5:400]*1e-6
                    r_0_vec
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec5 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for T_s = 18:32
                        T_s
                        calcCk_1drop;
                        Ck_vec5(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat5','Ck_vec5');
                
                T_s = 20;
            catch
                continue
            end
        case 6
            Andreas05_Fig11_in
            try
                % ------------- vary S only  ------------------------------------------------------
                count =1;
                
                for S =linspace(30,40,11)
                    S
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec6 maxIt maxEr_req maxEr_uf tau_r_exact
                    calcCk_1drop;
                    Ck_vec6(count)=C_k(1);
                    count = count+1;
                end
                save('Ck_mat6','Ck_vec6');
            catch
                continue
            end
        case 7
            Andreas05_Fig11_in
            try
                % ------------- vary S with U10  ------------------------------------------------------
                count =1;
                for S =linspace(30,40,11)
                    S
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec7 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for U10 = 10:5:80
                        calcCk_1drop;
                        Ck_vec7(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat7','Ck_vec7');
                U10 = 20;
            catch
                continue
            end
            
        case 8
            Andreas05_Fig11_in
            try
                % ------------- vary S with RH ------------------------------------------------------
                count =1;
                for S =linspace(30,40,11)
                    S
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec8 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for RH = 70:90
                        calcCk_1drop;
                        Ck_vec8(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat8','Ck_vec8');
                RH = 90;
            catch
                continue
            end
        case 9
            Andreas05_Fig11_in
            try
                % ------------- vary S with T_s  ------------------------------------------------------
                count =1;
                for S =linspace(30,40,11)
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec9 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for T_s = 18:32
                        calcCk_1drop;
                        Ck_vec9(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat9','Ck_vec9');
                T_s = 20;
                
                S = 34;%psu
            catch
                continue
            end
        case 10
            Andreas05_Fig11_in
            try
                
                % ------------- vary U10 only  ------------------------------------------------------
                count =1;
                
                for U10 = 10:5:80
                    U10
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec10 maxIt maxEr_req maxEr_uf tau_r_exact
                    calcCk_1drop;
                    Ck_vec10(count)=C_k(1);
                    count = count+1;
                end
                save('Ck_mat10','Ck_vec10');
                
            catch
                continue
            end
            
        case 11
            Andreas05_Fig11_in
            try
                % ------------- vary U10 with RH  ------------------------------------------------------
                count =1;
                
                for U10 = 10:5:80
                    U10
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec11 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for RH = 70:90
                        calcCk_1drop;
                        Ck_vec11(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat11','Ck_vec11');
                RH = 90;
            catch
                continue
            end
        case 12
            Andreas05_Fig11_in
            try
                % ------------- vary U10 with T_s  ------------------------------------------------------
                count =1;
                
                for U10 = 10:5:80
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec12 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for T_s = 18:32
                        calcCk_1drop;
                        Ck_vec12(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                save('Ck_mat12','Ck_vec12');
                
                T_s =20;
                
                U10 = 20;
            catch
                continue
            end
            
        case 13
            Andreas05_Fig11_in
            try
                % ------------- vary RH  only  ------------------------------------------------------
                count =1;
                
                for RH = 70:90
                    RH
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec13 maxIt maxEr_req maxEr_uf tau_r_exact
                    calcCk_1drop;
                    Ck_vec13(count)=C_k(1);
                    count = count+1;
                end
                save('Ck_mat13','Ck_vec13');
                
                
                count =1;
            catch
                continue
            end
        case 14
            Andreas05_Fig11_in
            try
                % ------------- vary RH  with T_s  ------------------------------------------------------
                
                for RH = 70:90
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec14 maxIt maxEr_req maxEr_uf tau_r_exact
                    countVar = 1;
                    for T_s = 18:32
                        calcCk_1drop;
                        Ck_vec14(count,countVar)=C_k(1);
                        countVar=countVar+1;
                    end
                    count = count+1;
                end
                
                save('Ck_mat14','Ck_vec14');
                RH = 90;
            catch
                continue
            end
        case 15
            Andreas05_Fig11_in
            try
                % ------------- vary  T_s only ------------------------------------------------------
                
                
                count =1;
                
                for T_s = 18:32
                    clearvars -except count r_0_vec  RH T_a T_s S p0 U10 Ck_vec15 maxIt maxEr_req maxEr_uf tau_r_exact
                    calcCk_1drop;
                    Ck_vec15(count)=C_k(1);
                    count = count+1;
                end
                save('Ck_mat15','Ck_vec15');
            catch
                continue
            end
    end
    
end
