function T = Stratified_T (T_inf, y, grad)
% calculate the temperature in a statified ambience
% T_inf is the reference temperature
% grad is the temperature gradient
% y is the vertical direction

T = T_inf + grad * y;

end