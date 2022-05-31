function f = polyfunction(p, x)
% calculate the function of the polyfit
% p is the coefficient of the function and x is the variable

L = length(p);
f = 0;
for i = 1: L
    f = f + p(i) * x^(L - i);    % iteration
end

end