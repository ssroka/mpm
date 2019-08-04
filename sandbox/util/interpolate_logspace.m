function [interpolant] = interpolate_logspace(x_vec,y_vec,x_val)


i_star = find(x_vec<x_val,1,'last');%index of the lower bound
slope  = (log10(y_vec(i_star+1))-log10(y_vec(i_star)))/(log10(x_vec(i_star+1))-log10(x_vec(i_star)));
interpolant  = 10.^(log10(y_vec(i_star))+slope*(log10(x_val)-log10(x_vec(i_star))));







end