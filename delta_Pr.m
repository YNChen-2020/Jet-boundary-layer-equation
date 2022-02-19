function delta_0 = delta_Pr(Pr)
% calculate the initial delta based on Pr

% the data from the table 1 in gersten-1980 paper
Pr_fit = [0.01; 0.1; 0.7; 1; 7; 10; 100;];
delta_fit = [39.370; 4.7192; 1.2241; 1.0; 0.36; 0.3004; 0.0945];

delta_0 = interp1(Pr_fit, delta_fit, Pr, 'linear');   % use the linear interpolation function

end
