function Rou = Rou_T(T)
% calculate the density related to the temperature (Tanaka - 2001)

a(1) = -3.983035;
a(2) = 301.797;
a(3) = 522528.9;
a(4) = 69.349;
a(5) = 999.975;
Rou_fun = @(t) a(5) * (1 - (t + a(1))^2 * (t + a(2)) / (a(3) * (t + a(4))));    % create the function handle
Rou = arrayfun(Rou_fun, T);

end
