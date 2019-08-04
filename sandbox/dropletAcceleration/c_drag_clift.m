function [c_d] = c_drag_clift(u,d,nu)

Re = u*d/nu;
w = log10(Re);

if Re < 0.01
    c_d = 3/16+24/Re;
elseif Re>0.01 && Re<=20
c_d = 24/Re*(1+0.1315*Re^(0.82-0.05*w));
elseif Re>20 && Re <260
    c_d = 24/Re*(1+0.1935*Re^(0.6305));
elseif Re >= 260 && Re <1500
    c_d = 10^(1.6435-1.1242*w+0.1558*w^2);
elseif Re >= 1500 && Re < 12000
    c_d = 10^(-2.4571+2.5558*w-0.929*w^2+0.1049*w^3);
elseif Re >= 12000 && Re < 44000
    c_d = 10^(-1.9181+0.6370*w-0.0636*w^2);
else
    Re
    error('Code more please')
end





















end







