%% load the data of experiment and simulation
clear; close all

% load the data mat file (case 12)
number_case = 12;
cd 'C:\Users\DELL\OneDrive\boundary-layer equations\simulation_data\case12\angle';
namelist = dir ('**\*.mat');
cd 'C:\Users\DELL\OneDrive\boundary-layer equations\code';
len = length(namelist);
namelist ((1 : len), :) = namelist([2:(len-1)/2 1 len (len+1)/2 len-1:-1:(len+3)/2], :);

% load the data of the simulation (different angles in case 12)
[s_s, x_s, y_s, alpha_s, U_s, width_s] = deal(cell(len, 1));
for i = 1 : len
    file_folder = namelist(i).folder;
    file_name = namelist(i).name;
    file_path = fullfile(file_folder, file_name);
    rawdata_s = struct2cell(load(file_path, 'S_distance',...
        'R_nond', 'R_d', 'x_real', 'y_real', 'U_initial'));
    % 'R_nond' and 'R_d' are cell variables and contain five differnet single variables
    % R_nond contains non-dimentional parameters Um, bu, delta, thetam, alpha
    % R_d contains dimentional parameters U, T, R_u, R_t, alpha
    
    data_s.S_distance = rawdata_s{1};    % create the structure matrix
    data_s.R_nond = rawdata_s{2};
    data_s.R_d = rawdata_s{3};
    data_s.x_real = rawdata_s{4};
    data_s.y_real = rawdata_s{5};
    data_s.U_initial = rawdata_s{6};
    
    % the non-dimensional parameters
    R_nond.Um = rawdata_s{2}(:,1);
    R_nond.bu = rawdata_s{2}(:,2);
    R_nond.delta = rawdata_s{2}(:,3);
    R_nond.thetam = rawdata_s{2}(:,4);
    R_nond.alpha = rawdata_s{2}(:,5);
    
    % the dimensional parameters
    R_d.U = rawdata_s{3}(:,1);
    R_d.T = rawdata_s{3}(:,2);
    R_d.R_u = rawdata_s{3}(:,3);
    R_d.R_t = rawdata_s{3}(:,4);
    R_d.alpha = rawdata_s{3}(:,5);
    
    % summarize the data
    s_s{i} = data_s.S_distance;
    x_s{i} = data_s.x_real;
    y_s{i} = data_s.y_real;
    alpha_s{i} = R_nond.alpha;
    U_s{i} = R_d.U;
    width_s{i} = R_d.R_u;

end

% the characteristic length scales
xls_filename = 'C:\Users\DELL\OneDrive\boundary-layer equations\工况整理.xlsx';
L_BMN = xlsread(xls_filename, 'N20:N31') ./ 100;     % change to the SI units: m

%% compare the jet centerline between simulation and experiment
color = [{[0 0.4470 0.7410]};{[0.8500 0.3250 0.0980]};{[0.9290 0.6940 0.1250]};...
    {[0.4940 0.1840 0.5560]};{[0.4660 0.6740 0.1880]};{[0.3010 0.7450 0.9330]};...
    {[0.6350 0.0780 0.1840]};{'r'};{'m'}];
[diff_ys, q] = deal(cell(len, 1));
length_c = L_BMN(number_case);

% plot the positive angle result
figure('units','normalized','position',[0.2 0.2 0.55 0.7]);
for i = 1 : (len-1)/2
    diff_ys{i} = diff(y_s{i});         % calcualte the diff trajectory in the simulation
    q{i} = plot(x_s{i} ./ length_c, y_s{i} ./ length_c, '--', 'color',...
        color{i}, 'lineWidth', 3.5);
    hold on
end

% plot the angle zero result
i = (len-1)/2 + 1;
diff_ys{i} = diff(y_s{i});
q{i} = plot(x_s{i} ./ length_c, y_s{i} ./ length_c, 'k-.', 'lineWidth', 3.5);

% plot the negative angle result
for i = 1+(len+1)/2 : len
    diff_ys{i} = diff(y_s{i});         % calcualte the diff trajectory in the simulation
    q{i} = plot(x_s{i} ./ length_c, y_s{i} ./ length_c, 'color',...
        color{end + 1 - (i - (len+1)/2)}, 'lineWidth', 3.5);
end

% set the figure parameters
set(gca, 'FontSize', 20);
legend([{'-17\pi/36'},{'-4\pi/9'},{'-5\pi/12'},{'-7\pi/18'},{'-13\pi/36'},{'-\pi/3'},...
    {'-\pi/4'},{'-\pi/6'},{'-\pi/12'},{'0'},{'\pi/12'},{'\pi/6'},{'\pi/4'},{'\pi/3'},...
    {'13\pi/36'},{'7\pi/18'},{'5\pi/12'},{'4\pi/9'},{'17\pi/36'}],...
    'NumColumns', 6, 'location', 'northeast', 'FontSize', 18);
xlim([0 2.5]);ylim([-2.3 3.3])
xlabel('$x^{*}/l^{*}_{12}$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y^{*}/l^{*}_{12}$', 'Interpreter', 'latex', 'FontSize', 24);
filepath = 'C:\Users\DELL\OneDrive\boundary-layer equations\verification\jet_angle_trajectory';
% export_figimage(strcat(filepath, '\angle_trajectory'), 'fig', '-jpg', 4);

%% calculate the max centerline trajectory and neutural centerline trajectory
[max_centerline_s, neu_centerline_s, max_point] = deal(zeros(len, 1));
L = length(diff_ys{1});
[point, k, b, turning_point]= deal(cell(len, 1));

% calculate the position of turning points
for i = 1 : len
    m = 1;
    for j = 1 : L-1
        a1 = diff_ys{i}(j);
        a2 = diff_ys{i}(j + 1);
        multi = a1 * a2;       % calculate the product of the number a1 and a2
        if multi < 0
            point{i}(m, 1) = j;
            m = m + 1;
        end
    end
end

% remove the wrong points and restore the real points
for i = 1 : len
    b{i} = diff(point{i});
    turning_point{i}(1, 1) = 1;
    n = 1;
    for j = 1 : length(b{i})
        if b{i}(j) > 2
            turning_point{i}(n+1, 1) = j + 1;      % store the value of the position
            n = n + 1;        
        end
    end
end

% calculate the max centerline and neutural centerline height
for i = 1 : len
     turn_point_1 = point{i}(turning_point{i}(1));          % find the first turning point
     turn_point_2 = point{i}(turning_point{i}(2));          % find the second turning point
     turn_point_3 = point{i}(turning_point{i}(3));          % find the third turning point
     
     [max_centerline_s(i), max_point(i)]= max(y_s{i}(turn_point_1 : turn_point_2));    % find the maximum height
     neu_centerline_s(i) = mean(y_s{i}(turn_point_2 : turn_point_3));    % find the neutral height
end

max_centerline_s(max_centerline_s < 0) = 0;     % set the negative maximum centerline height to be zero
max_centerline_s = max_centerline_s ./ length_c;
neu_centerline_s = neu_centerline_s ./ length_c;

%% plot the maximum height realted to the angle
% fit the maximum centerline height and the neutral height (with tanh(x))
x = [-pi*17/36; -pi*4/9; -pi*5/12; -pi*7/18; -pi*13/36; -pi/3; -pi/4; -pi/6; -pi/12; 0;...
    pi/12; pi/6; pi/4; pi/3; pi*13/36; pi*7/18; pi*5/12; pi*4/9; pi*17/36];
ft = fittype('a * tanh(x + b)', 'independent', 'x', 'coefficients', {'a', 'b'});
% ft = fittype('- a / (x + b) + c ', 'independent', 'x', 'coefficients', {'a', 'b', 'c'});
fo = fitoptions('Method', 'NonlinearLeastSquares');
[f_max, gof_max]= fit(x, max_centerline_s, ft, fo);       % fit the maximum height
[f_neu, gof_neu]= fit(x, neu_centerline_s, ft, fo);       % fit the neutral height
x_fit = (-pi / 2 : 0.1 : x(end))';
y_max = f_max(x_fit);
y_neu = f_neu(x_fit);

% plot the figure
figure('units','normalized','position',[0.2 0.2 0.4 0.6]);
plot(x, max_centerline_s, '*', 'MarkerSize', 18, 'LineWidth', 2.5); hold on
plot(x, neu_centerline_s, 'o', 'MarkerSize', 18, 'LineWidth', 2.5);
plot(x_fit, y_max, 'LineStyle', '-', 'linewidth', 4);
plot(x_fit, y_neu, 'LineStyle', '-', 'linewidth', 4);

xlim([-pi*3/5 pi/2]);ylim([-0.175 0.22] ./ length_c);
xticks([-pi*3/4, -pi/2, -pi/4, 0, pi/4, pi/2]);
xticklabels({'-3\pi/4','-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2'});
set(gca, 'FontSize', 22);
xlabel('$angle$','Interpreter','latex','FontSize',24);
ylabel('$y^{*}/l^{*}_{12}$','Interpreter','latex','FontSize',24);
legend('$H_{c-m}$', '$H_{n}$', 'Function 1', 'Function 2', 'Position',...
    [0.165, 0.786, 0.364, 0.117], 'FontSize', 18, 'Interpreter', 'latex',...
    'Orientation', 'vertical', 'NumColumns', 2);
% export_figimage(strcat(filepath, '\angle_height'), 'fig', '-jpg', 4);

%% plot the ratio of the height with regard to the result of angle zero
% fit the ratio function for maximum height and neutral height
max_norm = max_centerline_s ./ max_centerline_s((len + 1) / 2);
neu_norm = neu_centerline_s ./ neu_centerline_s((len + 1) / 2);
[f_max_norm, gof_max_norm]= fit(x, max_norm, ft, fo);
[f_neu_norm, gof_neu_norm]= fit(x, neu_norm, ft, fo);
y_max_norm = f_max_norm(x_fit);
y_neu_norm = f_neu_norm(x_fit);

% find the average fit function
f_norm = f_max_norm;
f_norm.a = (f_max_norm.a + f_neu_norm.a) / 2;
f_norm.b = (f_max_norm.b + f_neu_norm.b) / 2;
% f_norm.c = (f_max_norm.c + f_neu_norm.c) / 2;
y_norm = f_norm(x_fit);

figure('units','normalized','position',[0.2 0.2 0.4 0.6]);
plot(x, max_norm, '*', 'MarkerSize', 18, 'LineWidth', 2.5); hold on
plot(x, neu_norm, 'o', 'MarkerSize', 18, 'LineWidth', 2.5);
% plot(x_fit, y_neu_norm, 'LineStyle', '-', 'linewidth', 4);
% plot(x_fit, y_max_norm, 'LineStyle', '-', 'linewidth', 4);
plot(x_fit, y_norm, 'k', 'LineStyle', '-', 'linewidth', 4);

% set the figure parameters
xlim([-pi*3/5 pi/2]);ylim([-1.5 1.75]);
xticks([-pi*3/4, -pi/2, -pi/4, 0, pi/4, pi/2]);
xticklabels({'-3\pi/4','-\pi/2', '-\pi/4', '0', '\pi/4', '\pi/2'});
set(gca, 'FontSize', 22);
xlabel('$angle$','Interpreter','latex','FontSize',24);
ylabel('$R(\alpha)$','Interpreter','latex','FontSize',24);
legend('$H_{c-m}/H_{c-m}(0)$', '$H_{n}/H_{n}(0)$', '$R(\alpha)$', 'Position', ...
    [0.158, 0.716, 0.273, 0.181], 'FontSize', 20, 'Interpreter', 'latex',...
    'Orientation', 'vertical', 'NumColumns', 1);
% export_figimage(strcat(filepath, '\angle_function'), 'fig', '-jpg', 4);





