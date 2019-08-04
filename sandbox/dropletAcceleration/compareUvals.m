
%{
load /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/dropletAcceleration/plot2004.mat u_vec acel_time_vec
u_vec1 = u_vec;
acel_time_vec1 = acel_time_vec;

load /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/dropletAcceleration/plot2004_2.mat u_vec acel_time_vec
u_vec2 = u_vec;
acel_time_vec2 = acel_time_vec;

load /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/Calc_CK_vars.mat u_vec acel_time_vec
u_vec3 = u_vec;
acel_time_vec3 = acel_time_vec;

plot(acel_time_vec1,u_vec1,'r')
hold on
plot(acel_time_vec2,u_vec2,'k--')
plot(acel_time_vec3,u_vec3,'b')
%}


load /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/dropletAcceleration/plot204_jbode.mat  ic
ic

load /Users/ssroka/Documents/MATLAB/ASF/calcCorff_1drop/sandbox/c_ck.mat  ic
ic



























