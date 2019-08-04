%{
This script defines constants useful for the Andreas 2005 microphysical
model.

% 1: Defines useful anonymous functions
    Celsius2Kelvin, Kelvin2Celsius, Pa2mb, mb2Pa

% 2: Universal Constants
    g, R

% 3: Constants specific to air
    c_p, c_pd, M_a,


%}


% a few useful anonymous functions
helpingAnonFxns;


%% -------------------------------------------------------
%           Universal Constants
% --------------------------------------------------------
g = 9.81;    % [m s^-2] gravity
R = 8.31447; % [J mol^-1 K^-1] universal gas constant : see Andreas 2005 Symbols section

%% -------------------------------------------------------
%           Air
% --------------------------------------------------------
c_pd = 1.006e3;  % [J kg^-1 K^-1] specific heat of dry air Andreas 1995 section 2
M_a = 0.0289644; % [kg mol^-1]


%% -------------------------------------------------------
%           Salt
% --------------------------------------------------------

M_s = 0.058443; % [kg/mol] molecular weight of salt; Andreas 2005 symbols
nu  = 2;         % the number of ions NaCl dissociates into

%% -------------------------------------------------------
%            Water
% --------------------------------------------------------
M_w = 0.018015; % [kg/mol] molecular weight of water; Andreas 2005 Symbols










