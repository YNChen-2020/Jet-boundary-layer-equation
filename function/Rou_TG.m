function Rou_G = Rou_TG(T)
% calculate the gradient of density related to the temperature (Tanaka - 2001)

syms t
a(1) = -3.983035;
a(2) = 301.797;
a(3) = 522528.9;
a(4) = 69.349;
a(5) = 999.975;
Rou_fun = a(5) * (1 - (t + a(1))^2 * (t + a(2)) / (a(3) * (t + a(4))));    % create the density funcition
G_Rou_T = matlabFunction(diff(Rou_fun, t));    % create the gradient of the density function
Rou_G = G_Rou_T(T);

end
