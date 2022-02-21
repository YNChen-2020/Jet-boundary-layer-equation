function T_G = T_Rou(Rou)
% calculate the temperature based on the density (Tanaka - 2001)

syms t
a(1) = -3.983035;
a(2) = 301.797;
a(3) = 522528.9;
a(4) = 69.349;
a(5) = 999.975;
Rou_fun = a(5) * (1 - (t + a(1))^2 * (t + a(2)) / (a(3) * (t + a(4))));    % create the density funcition

fun = matlabFunction(Rou_fun - Rou);    % create the gradient of the density function
options = optimoptions('fsolve','Display','off');
T_G = fsolve (fun, 50, options);

end
