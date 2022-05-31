function  result = Rou_TS(T0, S0)
%% calculate the temperature and salinity based on the density (jirka - 2004)
syms T S
% calculate the density related to temperature (jirka - 2004)
Rou_ref = 999.842594;
T_coe(1) = 6.793952e-2;
T_coe(2) = -9.095290e-3;
T_coe(3) = 1.001685e-4;
T_coe(4) = -1.120083e-6;
T_coe(5) = 6.536332e-9;
Rou_T = Rou_ref + T_coe(1)*T + T_coe(2)*T^2 + T_coe(3)*T^3 ...
    + T_coe(4)*T^4 + T_coe(5)*T^5;    % the density function only related to temperature

% calculate the density related to temperature and salinity (jirka - 2004)
S_coe1_T(1) = 8.24493e-1;
S_coe1_T(2) = -4.0899e-3;
S_coe1_T(3) = 7.6438e-5;
S_coe1_T(4) = -8.2467e-7;
S_coe1_T(5) = 5.3875e-9;
S_coe(1) = S_coe1_T(1) + S_coe1_T(2)*T + S_coe1_T(3)*T^2 ...
    + S_coe1_T(4)*T^3 + S_coe1_T(5)*T^4;         % first salinity coefficient

S_coe2_T(1) = -5.72466e-3;
S_coe2_T(2) = 1.0227e-4;
S_coe2_T(3) = -1.6546e-6;
S_coe(2) = S_coe2_T(1) + S_coe2_T(2)*T + S_coe2_T(3)*T^2;    % second salinity coefficient

S_coe(3) = 4.8314e-4;       % third salinity coefficient

Rou_TS = Rou_T + S_coe(1)*S + S_coe(2)*S^(3/2)...
    + S_coe(3)*S^2;    % the density function related to temperature and density

Rou_fun = matlabFunction(Rou_TS, 'Vars', {T, S});
result = Rou_fun(T0, S0);


end





