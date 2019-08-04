    % a few useful anonymous functions
    Celsius2Kelvin =@(degC)degC+273.15;
    Kelvin2Celsius =@(K)K-273.15;
    Pa2mb =@(Pa)0.01*Pa;
    mb2Pa =@(mb)100*mb;