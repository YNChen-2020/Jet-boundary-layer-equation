function Mu = Mu_T(T)
% calculate the dynamic viscosity related to the temperature (Lide - 2015)

b(1) = 0.02939 * 10^-3;             % parameter one, SI units: Pa
b(2) = 507.88;                      % parameter two, SI units: K
b(3) = 149.3;                       % parameter three, SI units: K
Mu_T = @(t) b(1) * exp(b(2) / (t + 273.15 - b(3)));      % temperature should add 273.15K
Mu = Mu_T(T);

end