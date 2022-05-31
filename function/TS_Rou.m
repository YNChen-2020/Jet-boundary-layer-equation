function  result = TS_Rou(varargin)
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

%% calculate the temperature (gradient) based on salinity or salinity (gradient) based on temperature
Rou = varargin{1};
Rou_fun = matlabFunction(Rou_TS - Rou, 'Vars', {T, S});
% options = optimoptions('solve','Display','off');
dimension = varargin{3};

if nargin == 3
    % no gradient
    if dimension == 'T'              % judeg the type of input parameters
        T0 = varargin{2};
        equ_S = Rou_fun(T0, S);
        result_S = double(solve(equ_S, S, 'Real', true));    % calculate the salinity
        % delete the complex result in salinity
        index = real(result_S) > 0 & real(result_S) < 50;
        result = result_S(index);
    elseif dimension == 'S'          % judeg the type of input parameters
        S0 = varargin{2};
        equ_T = Rou_fun(T, S0);
        result_T = double(solve(equ_T, T, 'Real', true));    % calculate the temperature
        % delete the complex result in tempreature
        index = real(result_T) > 0 & real(result_T) < 100;
        result = result_T(index);
    end
    
    % gradient
elseif nargin == 5
    if dimension == 'T'              % judeg the type of input parameters
        T0 = varargin{2};
        S_tan = varargin{5};
        equ_S = Rou_fun(T0, S);
        equ_S_grad = matlabFunction(diff(equ_S, S));    % calculate the gradient function
        result = equ_S_grad(S_tan);
    elseif dimension == 'S'          % judeg the type of input parameters
        S0 = varargin{2};
        T_tan = varargin{5};
        equ_T = Rou_fun(T, S0);
        equ_T_grad = matlabFunction(diff(equ_T, T));    % calculate the gradient function
        result = equ_T_grad(T_tan);
    end
end


end













