% plot sensitivity study
clear;clc;close all

addpath ~/Documents/MATLAB/util/othercolor/
colormap(othercolor('RdYIGn4'))

r_0_vec = [10:5:400]*1e-6;
S =linspace(30,40,11);
U10 = 10:5:80;
RH  = 70:90;
T_s= 18:32;

%% ------------- vary r_0 [m] -------------------------------------------------

load Ck_mat1
figure(1)
subplot(2,3,1)
plot(r_0_vec,Ck_vec1)
xlabel('r_0 [m]')
ylabel('C_K')
sensitivity_study_plot_params

%% ------------- vary r_0 [m] with S ------------------------------------------
load Ck_mat2
figure(2)
subplot(3,3,1)
surf(r_0_vec,S,Ck_vec2')
xlabel('r_0 [m]')
ylabel('S [psu]')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary r_0 [m] with U10 ----------------------------------------
load Ck_mat3
figure(2)
subplot(3,3,2)
surf(r_0_vec,U10, Ck_vec3')
xlabel('r_0 [m]')
ylabel('U_{10} [m/s]')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary r_0 [m] with RH % -----------------------------------------
load Ck_mat4
figure(2)
subplot(3,3,3)
surf(r_0_vec,RH, Ck_vec4')
xlabel('r_0 [m]')
ylabel('RH %')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary r_0 [m] with T_s [deg C] ----------------------------------------
% load Ck_mat5
% figure(5)
% surf(r_0_vec,T_s [deg C], Ck_vec5')
% xlabel('r_0 [m]')
% ylabel('T_s [deg C]')
% zlabel('C_K')
% sensitivity_study_plot_params

%% ------------- vary S -------------------------------------------------
load Ck_mat6
figure(1)
subplot(2,3,2)
plot(S,Ck_vec6)
xlabel('S [psu]')
ylabel('C_K')
sensitivity_study_plot_params

%% ------------- vary S with U10 ----------------------------------------
load Ck_mat7
figure(2)
subplot(3,3,4)
surf(S,U10, Ck_vec7')
xlabel('S [psu]')
ylabel('U_{10} [m/s]')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary S with RH % ----------------------------------------
load Ck_mat8
figure(2)
subplot(3,3,5)
surf(S,RH, Ck_vec8')
xlabel('S [psu]')
ylabel('RH %')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary S with T_s [deg C] ----------------------------------------
load Ck_mat9
figure(2)
subplot(3,3,6)
surf(S,T_s, Ck_vec9')
xlabel('S [psu]')
ylabel('T_s [deg C]')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary U10 -------------------------------------------------

load Ck_mat10
figure(1)
subplot(2,3,3)
plot(U10,Ck_vec10)
xlabel('U_{10} [m/s]')
ylabel('C_K')
sensitivity_study_plot_params

%% ------------- vary U10 with RH % ----------------------------------------
load Ck_mat11
figure(2)
subplot(3,3,7)
surf(U10,RH, Ck_vec11')
xlabel('U_{10} [m/s]')
ylabel('RH %')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary U10 with T_s [deg C] ----------------------------------------
load Ck_mat12
figure(2)
subplot(3,3,8)
surf(U10,T_s, Ck_vec12')
xlabel('U_{10} [m/s]')
ylabel('T_s [deg C]')
zlabel('C_K')
sensitivity_study_plot_params


%% ------------- vary RH % only  ----------------------------------------
load Ck_mat13
figure(1)
subplot(2,3,4)
plot(RH,Ck_vec13)
xlabel('RH %')
ylabel('C_K')
sensitivity_study_plot_params

%% ------------- vary RH % with T_s [deg C] ----------------------------------------
load Ck_mat14
figure(2)
subplot(3,3,9)
surf(RH,T_s , Ck_vec14')
xlabel('RH %')
ylabel('T_s [deg C]')
zlabel('C_K')
sensitivity_study_plot_params

%% ------------- vary T_s [deg C] only  ----------------------------------------
load Ck_mat15
figure(1)
subplot(2,3,5)
plot(T_s ,Ck_vec15)
xlabel('T_s [deg C]')
ylabel('C_K')
sensitivity_study_plot_params











