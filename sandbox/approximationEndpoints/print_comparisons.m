
if exist('T_eq_exact','var') && print_flag
    fprintf('                         The value from Andreas 2005 is T_eq = %2.2f %sC\n',T_eq_exact,char(176))
    fprintf('     The full microphysical model value I computed is min(T) = %2.2f %sC\n',min_T,char(176));%char from ANSI char number
    fprintf('       The full microphysical model value I computed is T_eq = %2.2f %sC\n',T_eq_fmm,char(176));%char from ANSI char number
    fprintf('The estimated value (Andreas 2005 method) I computed is T_eq = %2.2f %sC\n\n',Kelvin2Celsius(T_eq),char(176));%char from ANSI char number
end

if exist('r_eq_exact','var') && print_flag
    fprintf('                         The value from Andreas 2005 is r_eq = %1.2f %sm\n',r_eq_exact*1e6,char(181))
    fprintf('       The full microphysical model value I computed is r_eq = %1.2f %sm\n',r_eq_fmm*1e6,char(181));%char from ANSI char number
    fprintf('The estimated value (Andreas 2005 method) I computed is r_eq = %1.2f %sm\n\n',r_eq*1e6,char(181));%char from ANSI char number
end

if exist('tau_T_exact','var')&& print_flag
    fprintf('                         The value from Andreas 2005 is tau_T = %1.3f s\n',tau_T_exact)
    fprintf('       The full microphysical model value I computed is tau_T = %1.3f s\n',tau_T_fmm);%char from ANSI char number
    fprintf('The estimated value (Andreas 2005 method) I computed is tau_T = %1.3f s\n\n',tau_T);%char from ANSI char number
end



if exist('tau_r_exact','var')&& print_flag
  fprintf('(I believe there is an error in)                        The value from Andreas 2005 is tau_r = %3.2f s\n',tau_r_exact)
    fprintf('(Since) The value extracted from the Andreas 2005 Fig 11 plot with WebPlotDigitizer is tau_r = %3.2f s\n',278.256)
    fprintf('                                      The full microphysical model value I computed is tau_r = %3.2f s\n',tau_r_fmm);%char from ANSI char number
    fprintf('                               The estimated value (Andreas 2005 method) I computed is tau_r = %3.2f s\n\n',tau_r);%char from ANSI char number
end