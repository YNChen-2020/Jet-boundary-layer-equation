%% load the data of experiment and simulation
clear; close all

% load the data mat file (case 12)
number_case = 12;
cd 'C:\Users\DELL\OneDrive\boundary-layer equations\simulation_data\case12\momentum';
namelist = dir ('**\*.mat');
cd 'C:\Users\DELL\OneDrive\boundary-layer equations\code';
len = length(namelist);
namelist ((1 : len), :) = namelist([1:len/2-1 len-1 len/2:len-2 len], :);

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

%% compare the jet centerline of different momentum
color = [{[0 0.4470 0.7410]};{[0.8500 0.3250 0.0980]};{[0.9290 0.6940 0.1250]};...
    {[0.4940 0.1840 0.5560]};{[0.4660 0.6740 0.1880]};{[0.3010 0.7450 0.9330]};...
    {[0.6350 0.0780 0.1840]};{'r'};{'m'};{'c'}];
[diff_ys, q] = deal(cell(len, 1));
length_c = L_BMN(number_case);

% plot the momentum
figure('units','normalized','position',[0.2 0.2 0.5 0.7]);
for i = 1 : len/2 - 1
    diff_ys{i} = diff(y_s{i});         % calcualte the diff trajectory in the simulation
    q{i} = plot(x_s{i} ./ length_c, y_s{i} ./ length_c, '-', 'color',...
        color{i}, 'lineWidth', 5);
    hold on
end

i = len/2;
diff_ys{i} = diff(y_s{i});
q{i} = plot(x_s{i} ./ length_c, y_s{i} ./ length_c, 'k-.', 'lineWidth', 5);

for i = len/2 + 1 : len
    diff_ys{i} = diff(y_s{i});         % calcualte the diff trajectory in the simulation
    q{i} = plot(x_s{i} ./ length_c, y_s{i} ./ length_c, '--', 'color',...
        color{i}, 'lineWidth', 5);
    hold on
end

% set the figure parameters
set(gca, 'FontSize', 24);
legend([{'$0.2U_{0}$'},{'$0.4U_{0}$'},{'$0.6U_{0}$'},{'$0.8U_{0}$'},{'$1.0U_{0}$'},...
    {'$1.2U_{0}$'},{'$1.4U_{0}$'},{'$1.6U_{0}$'},{'$1.8U_{0}$'},{'$2.0U_{0}$'}],...
    'NumColumns', 6, 'location', 'northeast', 'FontSize', 24, 'Interpreter', 'latex');
xlim([0 3.8]);ylim([0 2.2])
xlabel('$x^{*}/l^{*}_{12}$', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('$y^{*}/l^{*}_{12}$', 'Interpreter', 'latex', 'FontSize', 28);
filepath = 'C:\Users\DELL\OneDrive\boundary-layer equations\verification\jet_momentum_trajectory';
% export_figimage(strcat(filepath, '\momentum_trajectory'), 'fig', '-jpg', 4);

%% calculate the max centerline trajectory and neutural centerline trajectory
[max_centerline_s, neu_centerline_s, max_point] = deal(zeros(len, 1));
[point, k, b, turning_point]= deal(cell(len, 1));

% calculate the position of turning points
for i = 1 : len
    m = 1;
    L = length(diff_ys{i});
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
max_centerline_sn = max_centerline_s ./ length_c;
neu_centerline_sn = neu_centerline_s ./ length_c;
ratio = neu_centerline_sn ./ max_centerline_sn;

%% plot the maximum height realted to the angle
% fit the maximum centerline height and the neutral height (with tanh(x))
x = [0.2; 0.4; 0.6; 0.8; 1.0; 1.2; 1.4; 1.6; 1.8; 2];
ft = fittype('a * x + b', 'independent', 'x', 'coefficients', {'a', 'b'});
fo = fitoptions('Method', 'NonlinearLeastSquares');
[f_max, gof_max]= fit(x, max_centerline_sn, ft, fo);       % fit the maximum height
[f_neu, gof_neu]= fit(x, neu_centerline_sn, ft, fo);       % fit the neutral height
x_fit = (x(1) : 0.1 : x(end))';
y_max = f_max(x_fit);
y_neu = f_neu(x_fit);

% plot the figure
figure('units','normalized','position',[0.2 0.2 0.4 0.6]);
plot(x, max_centerline_sn, '*', 'MarkerSize', 22, 'LineWidth', 3); hold on
plot(x, neu_centerline_sn, 'o', 'MarkerSize', 22, 'LineWidth', 3);
plot(x_fit, y_max, 'LineStyle', '-', 'linewidth', 5);
plot(x_fit, y_neu, 'LineStyle', '-', 'linewidth', 5);

xlim([x(1) x(end)]);
ylim([0.4 1.1 * max(max_centerline_sn)]);
xticks([0.4 0.8 1.2 1.6 2])
set(gca, 'FontSize', 24);
xlabel('$U_{jet}/U_{0}$','Interpreter','latex','FontSize',26);
ylabel('$y^{*}/l^{*}_{12}$','Interpreter','latex','FontSize',26);
legend('$H_{c-m}$', '$H_{n}$', 'Function 1', 'Function 2', 'Position',...
    [0.495, 0.777, 0.376, 0.117], 'FontSize', 22, 'Interpreter', 'latex',...
    'Orientation', 'vertical', 'NumColumns', 2);

fprintf('中心线最大高度的方程为 %.3f *x + %.3f b \n', f_max.a, f_max.b)
fprintf('中性浮力层高度的方程为 %.3f *x + %.3f b \n', f_neu.a, f_neu.b)

% export_figimage(strcat(filepath, '\momentum_height'), 'fig', '-jpg', 4);

%% calculate the real height related to the length scale
L_new = length_c ./ x.^(3/4);
figure('units','normalized','position',[0.2 0.2 0.4 0.6]);
plot(L_new, max_centerline_s ./ L_new, '*', 'MarkerSize', 22, 'LineWidth', 3); hold on
plot(L_new, neu_centerline_s ./ L_new, 'o', 'MarkerSize', 22, 'LineWidth', 3);
set(gca, 'FontSize', 24);
xlabel('$l^{*} \ (m)$','Interpreter','latex','FontSize',30);
ylabel('$y^{*} \ (m)$','Interpreter','latex','FontSize',30);

% save(strcat(filepath, '\momentum'));

%% calculate the real height related to the length scale
% figure('units','normalized','position',[0.2 0.2 0.4 0.6]);
% plot(x, max_centerline_s, '*', 'MarkerSize', 22, 'LineWidth', 3); hold on
% plot(x, neu_centerline_s, 'o', 'MarkerSize', 22, 'LineWidth', 3);
% set(gca, 'FontSize', 24);
% xlabel('$U_{jet}/U_{0}$','Interpreter','latex','FontSize',30);
% ylabel('$y^{*} \ (m)$','Interpreter','latex','FontSize',30);





