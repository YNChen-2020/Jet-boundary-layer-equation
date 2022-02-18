%% set the physical parameters
close all
clear; clc
addpath('C:\Users\DELL\OneDrive\boundary-layer equations\code\function')

% the density of water related to the temperature (Tanaka - 2001)
Rou_0 = Rou_T(4);           % density (4℃), SI units: Kg/m^(3) 

% the dynamic viscosity of water related to the tempreature (Lide - 2015)
Mu_0 = Mu_T(40);            % dynamic viscosity (20℃, water), SI units: N*S/m^(-2)

% the thermal conductivity of water related to the temperature (the difference is not much == constant, wiki)
C_p = 4184;                 % specifc heat, SI units: J/(Kg*K)

% the tspecifc heat of water related to the temperature (the difference is not much == constant, wiki)
K = 0.6089;                 % thermal conductivity (26.85℃, water), SI units: J/(m*s*K)

% set the initial parameters
Pr = C_p * Mu_0 / K;        % Prandtl number, SI units: 1

%% set the initial condition of buoyant jets (my experiment case)
% near the inlet nozzle, the "top-hat" assumption is applied here, the

% load the parameters of jets in different cases
xls_filename = 'C:\Users\DELL\OneDrive\boundary-layer equations\工况整理.xlsx';
Rou_jet = xlsread(xls_filename, 5, 'D3:D14');
Rou_a = xlsread(xls_filename, 5, 'E3:E14');
U_jet = xlsread(xls_filename, 5, 'N3:N14');      % the value of the velocity
Density_gradient = xlsread(xls_filename, 5, 'V3:V14');

case_number = 12;

U_initial = U_jet(case_number) * 2;
% U_initial = 0.152;
T_initial = T_Rou(Rou_jet(case_number));
T_inf = T_Rou(Rou_a(case_number));

% test the negative buoyancy force
% T_inf = T_Rou(Rou_jet(case_number));
% T_initial = T_Rou(Rou_a(case_number));


% U_initial = 0.07;           % initial velocity, SI units: m/s
% T_initial = 40;             % initial temperature, SI units: ℃
% T_inf = 4;                  % ambient temperature, SI units: ℃
D = 3.7 * 10^(-3);          % the diameter of the nozzle, SI units: m
T_grad = Density_gradient(case_number) / 0.2 * 100;         % temperature gradient, SI units: ℃

% D_set = [16.8; 0.78; 1.19];
% D = D_set(case_number) * 10^(-3);          % the diameter of the nozzle, SI units: m
g = 9.8;                    % gravity acceleration, SI units: m/s^2

% the reference characteristic parameters for the first iteration
% S_linspace = 0.1;
% S_limup = (0.01 : S_linspace : 5)';
S_linspace = 0.02;         % the final calculation
S_limup = (0.01 : S_linspace : 20)';
len = length(S_limup);
S_itn_lin = 11;              % the linspace in the iteration is 11

% the initial values for the reference parameters
[belta_itn, U0, delta_T, L] = deal(zeros(len + 1, 1));
E_T = abs(U_initial * (T_initial - T_inf) * D);        % thermal energy
K_j = U_initial^2 * D;                            % kinematic momentum
belta_itn(1) = - Rou_TG(T_initial) / Rou_0;       % the approximate slope of the density at temperature 40℃
U0(1) = (E_T * g * belta_itn(1))^(1/3);                           % characteristic velocity
delta_T(1) = (E_T * g * belta_itn(1))^(1/3) * (E_T / K_j);        % characteristic temperature difference
L(1) = (E_T * g * belta_itn(1))^(-2/3) * K_j;                     % characteristic length

% the initial values for the five non-dimensional variables
[Um_int, bu_int, delta_int, thetam_int, alpha_int, Ta_int]...
    = deal(zeros(len + 1,1));
Um_int(1) = U_initial / U0(1);
bu_int(1) = D / (2 * L(1));
delta_int(1) = 0.36;             % the temperature field range is about 30% of the velocity field (Table 1 - gersten-1980)
thetam_int(1) = (T_initial - T_inf) / delta_T(1);
thetam_int2(1) = (T_initial - T_inf) / delta_T(1);
alpha_int(1) = 0;                % in radians
Ta_int(1) = T_inf;

% the initial values for the four dimensional variables
[U, T, R_u, R_t] = deal(zeros(len + 1, 1));
U(1) = U_initial;
T(1) = T_initial;
R_u(1) = D / 2;
R_t(1) = R_u(1) * delta_int(1);

% pre-allocate the memory of other parameters
S_distance = zeros(len, 1);
S_distance(1) = 0;
y_dis = zeros(len + 1, 1);
y_dis(1) = 0;
[S, R] = deal(cell(len, 1));
[Ar, Pr_t] = deal(zeros(len, 1));

%% set the shape parameters 
% A1 to A8 with exact value (reslts from gersten - 1980, verfied by me)
% A = zeros(len + 1, 7);
A1 = @(delta_a) -0.77 * (tanh (1.05 / delta_a))^0.85 + 2.89;
A2 = 1.513;
A3 = 1.210;
A4 = @(delta_a) 2.12 * tanh (0.7 / delta_a)^0.7;
A5 = 0.940;
A6 = @(delta_a) -2.0 * (tanh (0.51 / delta_a))^0.9;
A7 = @(delta_a) 1.44 * (tanh (0.65 / delta_a))^1.4;
A8 = 2.269;

% A(1,1) = A1(delta_int(1));
% A(1,2) = A2;
% A(1,3) = A3;
% A(1,4) = A4(delta_int(1));
% A(1,5) = A5;
% A(1,6) = A6(delta_int(1));
% A(1,7) = A7(delta_int(1));
% A(1,8) = A8;

syms U_m(s) theta_m(s) b_u(s) delta(s) alpha(s) s Y

[Um_middle, bu_middle, delta_middle, thetam_middle, alpha_middle]...
    = deal(zeros(len,1));

for i = 1: len
%% solve the differential equations

% set reference characteristics
U0_itn = U0(i);
L_int = L(i);
delta_T_int = delta_T(i);

% resluts from Schlichting's book (Page 661) 
Pr_t = 0.84;                                      % turbulent Prandtl number
Nu_t = 0.036 * U0_itn * L_int * U_m * b_u;        % kinematic eddy viscosity

% % resluts from Gersten-1980
% Ar(i) = g * belta_itn(i) * E_T * (A8 / A2)^3 * Um_int(i)^(-3);
% co = 0.035 * (1 + 1.64 * Ar(i));
% Nu_t = co * U0_itn * L_int * U_m * b_u;
% Pr_t(i) = 0.5 * (1 + 1.02 * Ar(i));

% define the differential equations
D1 = diff(U_m^2 * b_u * A2, s);
D2 = diff(U_m^3 * b_u * A3 / 2, s);
D3 = diff(U_m * theta_m * b_u * delta * A4(delta), s);
D4 = diff(U_m^2 * theta_m * b_u^2 * delta^2 * A7(delta), s);
D5 = diff(alpha, s);

eqn1 = D1 == theta_m * b_u * delta * A1(delta) * sin(alpha);             % in radians
eqn2 = D2 == U_m * theta_m * b_u * delta * A4(delta) * sin(alpha)...     % in radians
        - (Nu_t / (U0_itn * L_int)) * (U_m^2 * A5 / b_u);
eqn3 = D3 == 0;
eqn4 = D4 == - (Nu_t / (Pr_t * U0_itn * L_int)) * U_m * theta_m * A6(delta);
eqn5 = D5 == (theta_m * delta * A1(delta) * cos(alpha)) / (U_m^2 * A2);  % in radians

% the initial condition of this iteration
Um_0 = Um_int(i);
bu_0 = bu_int(i);
delta_0 = delta_int(i);
thetam_0 = thetam_int(i);
alpha_0 = alpha_int(i);

% solve the problem using the ode45 function (numerical solutions)
[VF, Sbs] = odeToVectorField(eqn1, eqn2, eqn3, eqn4, eqn5);    % reduce the order of differential equations
odefcn = matlabFunction(VF, 'Vars',{s, Y});                    % convert symbolic expression to function handle or file

S_limdown = S_limup(i) - S_linspace;
S_span = linspace(S_limdown, S_limup(i), S_itn_lin);           % the calculation range
Iconds = [Um_0, bu_0, delta_0, thetam_0, alpha_0];             % the initial conditions for variables in "Sbs"
[S{i}, R{i}] = ode45(odefcn, S_span, Iconds);                  % using the function "ode45" to solve
S_distance(i+1) = S{i}(end);

%% update the initial values of five variables for the next calculation region
Um_middle(i) = R{i}(end,1);
bu_middle(i) = R{i}(end,2);
delta_middle(i) = R{i}(end,3);
thetam_middle(i) = R{i}(end,4);
alpha_middle(i) = R{i}(end,5);

% Um_int(i+1) = R{i}(end,1);
% bu_int(i+1) = R{i}(end,2);
% delta_int(i+1) = R{i}(end,3);
% % thetam_int(i+1) = R{i}(end,4);
% alpha_int(i+1) = R{i}(end,5);

%% calculate the real variables, update the ambient temperature and reference length
U(i+1) = Um_middle(i) * U0_itn;                          % real velocity
T(i+1) = thetam_middle(i) * delta_T_int +  Ta_int(i);         % real temperature
R_u(i+1) = bu_middle(i) * L_int;                         % real radius related to velocity
R_t(i+1) = delta_middle(i) * R_u(i+1);                   % real radius related to temperature

% U(i+1) = Um_int(i+1) * U0_itn;                          % real velocity
% % T(i+1) = thetam_int(i+1) * delta_T_int + Ta_int(i);   % real temperature
% T(i+1) = thetam_int(i+1) * delta_T_int + T_inf;         % real temperature
% R_u(i+1) = bu_int(i+1) * L_int;                         % real radius related to velocity
% R_t(i+1) = delta_int(i+1) * R_u(i+1);                   % real radius related to temperature

% calculate the ambient temperature at real veritical distance
y_diff = cumsum(sin(R{i}(:,5))) .* S_linspace * L_int / (S_itn_lin - 1);   % the difference in vertical distance
y_dis(i+1) = y_dis(i) + y_diff(end);
Ta_int(i+1) = Stratified_T(Ta_int(i), y_diff(end), T_grad);                % use the costume function "Stratified_T" 

% calculate the reference velocity, temperaure, and length for the next calculation region
belta_itn(i+1) = - Rou_TG(T(i+1)) / Rou_0;                           % the approximate slope of the density at temperature T_jet
% belta_itn(i+1) = belta_itn(i);
% E_T(i+1) = A4(delta_int(i+1)) * Um_int(i+1) * (T(i+1) - Ta_int(i)) * R_t(i+1);
% K_j(i+1) = A2 * Um_int(i+1)^2 * R_u(i+1);
U0(i+1) = (E_T * g * belta_itn(i+1))^(1/3);                          % characteristic velocity
delta_T(i+1) = (E_T * g * belta_itn(i+1))^(1/3) * (E_T / K_j);       % characteristic temperature difference
L(i+1) = (E_T * g * belta_itn(i+1))^(-2/3) * K_j;                    % characteristic length

%% update the initial values of five variables for the next calculation region
thetam_int(i+1) = (T(i+1) - Ta_int(i+1)) / delta_T(i+1);
Um_int(i+1) = U(i+1) / U0(i+1);
bu_int(i+1) = R_u(i+1) / L(i+1);
delta_int(i+1) = R_u(i+1) / R_t(i+1);
alpha_int(i+1) = alpha_middle(i);

end

% summarize the final results
R_nond = [Um_int, bu_int, delta_int, thetam_int, alpha_int];      % non-dimensional
R_d = [U, T, R_u, R_t, alpha_int];                                % dimensional

%% plot the nondimentional parameters (loglog)
figure('units','normalized','position',[0.2 0.2 0.42 0.55]);
loglog(S_distance, R_nond, 'Linewidth',5); grid                   % change to the Log-log scale plot
set(gca,'FontSize',14);
legend('$U_{m}$','$b_{u}$','$\delta$','$\theta_{m}$','$\alpha$', 'Position',...
    [0.702,0.173,0.137,0.285],'FontSize',16,'Interpreter','latex')
title('Simulation of the boundary-layer equations','FontSize',20);
xlim([S_limup(1) S_limup(end)])
ylim([min(min(R_nond)) S_limup(end)])
xlabel('S','FontSize',20);
ylabel('value','FontSize',20);
filepath = 'C:\Users\DELL\OneDrive\boundary-layer equations\simulation_data';
case_name = strcat('\case', num2str(case_number));
% export_figimage(strcat(filepath, case_name, '\non-dimensional parameters_loglog'),'fig','-jpg',4);

%% plot the dimentional parameters (loglog)
figure('units','normalized','position',[0.2 0.2 0.42 0.55]);
loglog(S_distance, R_d, 'Linewidth',5); grid                      % change to the Log-log scale plot
set(gca,'FontSize',18);
legend('$U$','$T$','$\delta_{u}$','$\delta_{t}$','$\alpha$', 'Position',...
    [0.183,0.581,0.137,0.285],'FontSize',20,'Interpreter','latex')
title('Simulation of the boundary-layer equations','FontSize',24);
xlim([S_limup(1) S_limup(end)])
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$value$', 'Interpreter', 'latex', 'FontSize', 24);
% export_figimage(strcat(filepath, case_name, '\dimensional parameters_loglog'),'fig','-jpg',4);

%% plot the jet trajectory
% calculate the distance in the horizontal distance and vertical distance
x_real = cumsum(cos(alpha_int)) .* S_linspace .* L;      % argument in radians
y_real = cumsum(sin(alpha_int)) .* S_linspace .* L;      % argument in radians
figure('units','normalized','position',[0.2 0.2 0.42 0.55]);
plot(x_real, y_real, 'Linewidth',3); grid
title('Jet trajectory','FontSize',24);
set(gca,'FontSize',18);
xlabel('$x$','FontSize',24,'Interpreter','latex');
ylabel('$y$','FontSize',24,'Interpreter','latex');
% export_figimage(strcat(filepath, case_name, '\jet trajectory'),'fig','-jpg',4);

%% plot the nondimentional parameters (Descartes)
figure('units','normalized','position',[0.2 0.2 0.4 0.55]);
plot(S_distance, R_nond, 'Linewidth',6); grid                   % change to the Log-log scale plot
set(gca,'FontSize',18);
legend('$U_{m}$','$b_{u}$','$\Delta$','$\theta_{m}$','$\alpha$', 'Position',...
    [0.262,0.596,0.137,0.285],'FontSize',20,'Interpreter','latex')
% xlim([S_limup(1) S_limup(end)])
% ylim([min(min(R_nond)) S_limup(end)])
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$value$', 'Interpreter', 'latex', 'FontSize', 24);
% export_figimage(strcat(filepath, case_name, '\non-dimensional parameters_Descartes'),'fig','-jpg',4);

%% plot the dimentional parameters (Descartes)
fig = figure('units','normalized','position',[0.2 0.2 0.4 0.55]);
no_T = [1 3 4 5];
color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];   %the color of the line
left_color = color(1, :);
right_color = color(4, :);
set(fig, 'defaultAxesColorOrder', [left_color; right_color]);
set(gca,'FontSize',18);

% plot the left axis
yyaxis left
p = plot(S_distance, R_d(:, no_T),'-' , 'Linewidth', 6 );
p(1).Color = color(1, :);
p(2).Color = color(2, :);
p(3).Color = color(3, :);p(3).LineStyle = ':';        % set the color and linestyle
p(4).Color = color(5, :);
ylim([-1 1])
ylabel('$value$','Interpreter','latex','FontSize',24);

% plot the right axis
yyaxis right
plot(S_distance, R_d(:, 2), 'Linewidth', 6, 'color', color(4, :)); grid
ylabel('$temperature$','Interpreter','latex','FontSize',24);
legend('$\bar u_{m}$','$\delta_{u}$','$\delta_{t}$','$\alpha$','$\bar T_{m}$', 'Position',...
    [0.159,0.845,0.500,0.0650],'FontSize',20,'Interpreter','latex',...
    'Orientation','horizontal')
xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 24);
% export_figimage(strcat(filepath, case_name, '\dimensional parameters_Descartes'),'fig','-jpg',4);

%% save the workspace
% save(strcat(filepath, case_name, '\', case_name));

% % write the case setting
% fileID = fopen(strcat(filepath, case_name, '\settings.txt'), 'w');
% fprintf(fileID, '计算域为%1.1f  %2.2f \n', 0, S_limup(end));
% fprintf(fileID, '计算间隔为%2.2f \n', S_linspace);
% fclose(fileID);
