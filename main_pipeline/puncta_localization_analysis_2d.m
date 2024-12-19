%% Analysis pipeline performed using puncta localization

% Sayan Seal, Jan 2024

% Add necessary paths; the following paths are relative to the current directory being the "main_pipeline/" subdirectory
addpath(genpath('../LLSM5DTools/'));
addpath(genpath('../util_functions/'));

%% Load data

% root directory containing the subdirectories for each video, which consists of the csv files for the detections and the corresponding masks
% For replicating the density scatter and frequency distribution plots in https://doi.org/10.1101/2023.12.22.573117,the dataset can be downloaded from
% https://doi.org/10.5061/dryad.w6m905qwm, and the rootdir path needs to be changed accordingly.
rootdir = '/home/sayanseal/doi_10_5061_dryad_w6m905qwm__v20240307/three_plane_stack_data/BaptaWT/140um/';

% directory to store the results
resultsdir = [rootdir 'results/'];
mkdir(resultsdir);

csv_list = dir(fullfile(rootdir, '**/*locs_filtered.csv')); % All csv files corresponding to specific filtering conditions obtained from performing puncta localization using ThunderSTORM 
save_str = 'All_WT_Bapta_140um_locs_all_frames'; % filename string for plots

lumen_list = dir(fullfile(rootdir, '**/*lumen*.tif')); % Lumen Edge Mask files
outer_list = dir(fullfile(rootdir, '**/*outer*.tif')); % Outer Edge Mask files
csv_filenames = {csv_list.name};
csv_filepaths = {csv_list.folder};
lumen_filenames = {lumen_list.name};
lumen_filepaths = {lumen_list.folder};
outer_filenames = {outer_list.name};
outer_filepaths = {outer_list.folder};

%% Extract necessary distance information

% Distances of puncta localizations from masks across all videos
lumen_dist_arr = []; 
outer_dist_arr = [];

per_exp_pair = {}; % cell array combining [Distances from Lumen Edge Mask, Distances from Outer Edge Mask] pairs for each video

count = 0; % total number of detections
num_det = []; % number of detections per video

% CDF of the distances from masks across all videos
lumen_cdf_vals = [];
outer_cdf_vals = [];

% Total number of pixels for each mask
lumen_count = 0;
outer_count = 0;

for f = 1:numel(csv_filenames)
    % Lumen Edge Mask
    lumen = parallelReadTiff([lumen_filepaths{f} filesep lumen_filenames{f}]); % read mask
    dtm_lumen = bwdist(logical(lumen)); % perform distance transform
    invdtm_lumen = bwdist(~logical(lumen)) * -1; % option to indicate pixels lying on the mask with negative distances; can be commented out.
    cmbdtm_lumen = dtm_lumen;
    cmbdtm_lumen(logical(lumen)) = invdtm_lumen(logical(lumen));
    cmbdtm_lumen(cmbdtm_lumen < 0) = 0; % distance of pixels lying on the mask set to 0

    % Outer Edge Mask
    outer = parallelReadTiff([outer_filepaths{f} filesep outer_filenames{f}]); % read mask
    dtm_outer = bwdist(logical(outer)); % perform distance transform
    invdtm_outer = bwdist(~logical(outer)) * -1; % option to indicate pixels lying on the mask with negative distances; can be commented out.
    cmbdtm_outer = dtm_outer;
    cmbdtm_outer(logical(outer)) = invdtm_outer(logical(outer));
    cmbdtm_outer(cmbdtm_outer < 0) = 0; % distance of pixels lying on the mask set to 0

    lumen_dist_arr_per_exp = [];
    outer_dist_arr_per_exp = [];

    lumen_count = lumen_count + nnz(lumen);
    outer_count = outer_count + nnz(outer);

    localization = readtable([csv_filepaths{f} filesep csv_filenames{f}], 'VariableNamingRule', 'preserve');
%     loc = localization(localization.frame == 1,:);  % uncomment if only the first frame is needed, and comment out the next line.
    loc = localization;

    % Convert to pixel space
    x = round(loc.("x [nm]")/108);
    y = round(loc.("y [nm]")/108);

    count = count + length(x);
    num_det = [num_det, length(x)];

    for n = 1:numel(x)
        dist_lumen = cmbdtm_lumen(y(n),x(n))*.108; % convert to microns
        lumen_dist_arr_per_exp = [lumen_dist_arr_per_exp, dist_lumen];
        lumen_dist_arr = [lumen_dist_arr, dist_lumen];

        dist_outer = cmbdtm_outer(y(n),x(n))*.108; % convert to microns
        outer_dist_arr_per_exp = [outer_dist_arr_per_exp, dist_outer];
        outer_dist_arr = [outer_dist_arr, dist_outer];

    end

    per_exp_pair = [per_exp_pair; [lumen_dist_arr_per_exp; outer_dist_arr_per_exp]];

    [lumen_f, lumen_x] = ecdf(lumen_dist_arr_per_exp);
    lumen_cdf_vals = [lumen_cdf_vals; interp1(lumen_x(2:end), lumen_f(2:end), 0:1:80)];

    [outer_f, outer_x] = ecdf(outer_dist_arr_per_exp);
    outer_cdf_vals = [outer_cdf_vals; interp1(outer_x(2:end), outer_f(2:end), 0:1:80)];
end

%% Store the CDF values per video of the distances from masks

% Correction for NaN values
for i=1:numel(lumen_cdf_vals(:, 1))
    lumen_k = find(~isnan(lumen_cdf_vals(i, 1:end)));
    for j=1:lumen_k(1)-1
        lumen_cdf_vals(i, j) = 0;
    end
    outer_k = find(~isnan(outer_cdf_vals(i, 1:end)));
    for j=1:outer_k(1)-1
        outer_cdf_vals(i, j) = 0;
    end
end
lumen_cdf_vals(isnan(lumen_cdf_vals)) = 1;
outer_cdf_vals(isnan(outer_cdf_vals)) = 1;

% Organize the cells
lumen_cdf_vals_final = [csv_filenames' num2cell(lumen_cdf_vals)];
outer_cdf_vals_final = [csv_filenames' num2cell(outer_cdf_vals)];

headers = {'Filename'};
for n=2:length(lumen_cdf_vals_final)
    headers{1,n} = [num2str(n-2) ' ' char(956) 'm'];
end

lumen_cdf_vals_final = [headers; lumen_cdf_vals_final];
outer_cdf_vals_final = [headers; outer_cdf_vals_final];

cdf_filename = 'All_WT_Bapta_140um_locs_all_frames_per_file.xlsx';
writecell(lumen_cdf_vals_final, [resultsdir cdf_filename], 'Sheet', 'Lumen');
writecell(outer_cdf_vals_final, [resultsdir cdf_filename], 'Sheet', 'Outer Edge');

%% Store the combined CDF values across all videos of the distances from masks

[lumen_comb_f, lumen_comb_x] = ecdf(lumen_dist_arr);
lumen_cdf_comb_vals = [interp1(lumen_comb_x(2:end), lumen_comb_f(2:end), 0:1:80)];

[outer_comb_f, outer_comb_x] = ecdf(outer_dist_arr);
outer_cdf_comb_vals = [interp1(outer_comb_x(2:end), outer_comb_f(2:end), 0:1:80)];

% Correction for NaN values
lumen_comb_k = find(~isnan(lumen_cdf_comb_vals(1:end)));
for j=1:lumen_comb_k(1)-1
    lumen_cdf_comb_vals(j) = 0;
end
outer_comb_k = find(~isnan(outer_cdf_comb_vals(1:end)));
for j=1:outer_comb_k(1)-1
    outer_cdf_comb_vals(j) = 0;
end
lumen_cdf_comb_vals(isnan(lumen_cdf_comb_vals)) = 1;
outer_cdf_comb_vals(isnan(outer_cdf_comb_vals)) = 1;

% Organize the cells
headers = {};
for n=1:length(lumen_cdf_comb_vals)
    headers{1,n} = [num2str(n-1) ' ' char(956) 'm'];
end

lumen_cdf_comb_vals_final = [headers; num2cell(lumen_cdf_comb_vals)];
outer_cdf_comb_vals_final = [headers; num2cell(outer_cdf_comb_vals)];

cdf_comb_filename = 'All_WT_Bapta_140um_locs_all_frames_combined.xlsx';
writecell(lumen_cdf_comb_vals_final, [resultsdir cdf_comb_filename], 'Sheet', 'Lumen');
writecell(outer_cdf_comb_vals_final, [resultsdir cdf_comb_filename], 'Sheet', 'Outer Edge');

%% Density Scatter Plot of distances from Lumen Edge and Outer Edge Masks

fig = figure(); 
fig.Position = [1000, 1000, 1200, 1200];
densityScatterChart(double(lumen_dist_arr), double(outer_dist_arr), 'Colormap', jet); % marker size can be adjusted by modifying the 'SizeData' parameter in "../util_functions/densityScatterChart.m"
xlim([-1, 80]);
ylim([-1, 80]);
xticks(0:10:80);
yticks(0:10:80);
clim([0 2e-3]);
xlabel('Distance from Lumen Edge Mask (\mum)');
ylabel('Distance from Outer Edge Mask (\mum)');
title('Active PIEZO1-HaloTag Localizations');
% set(gca, 'XAxisLocation', 'top');
set(gca, 'PlotBoxAspectRatio', [1,1,1]);
set(gca, 'FontName', 'Arial', 'FontSize', 23);
set(gca, 'Position', [.12,.15,.75,.75]);
save_filename = ['Density_Scatter_' save_str  '.svg'];
saveas(gcf, [resultsdir save_filename]);

%% Frequency Distribution Plots of distances from Lumen Edge

fig = figure(); 
fig.Position = [1000, 1000, 1200, 1200];
GU_stairsHistPlot({lumen_dist_arr}, 'BinSize', 2, 'ShowECDF', true, 'LineWidth', 2,  'xlabel', 'Distance from Lumen (\mum)', 'ylabel', 'Relative Frequency', 'FontSize', 23, 'CreateFigure', false); % color can be adjusted by modifying the 'FaceColor', 'EdgeColor' and 'Color' parameters in "../util_functions/GU_stairsHistPlot.m"
xlim([-1, 80]);
xticks(0:10:80);
pbaspect([1,1,1]);
set(gca, 'Position', [.14,.12,.75,.75]);
pbaspect([1,1,1]);
box on
save_filename = ['Density_Lumen_' save_str  '.svg'];
saveas(gcf, [resultsdir save_filename]);

%% Frequency Distribution Plots of distances from Outer Edge

fig = figure(); 
fig.Position = [1000, 1000, 1200, 1200];
GU_stairsHistPlot({outer_dist_arr}, 'BinSize', 2, 'ShowECDF', true, 'LineWidth', 2,  'xlabel', 'Distance from Outer Edge (\mum)', 'ylabel', 'Relative Frequency', 'FontSize', 23, 'CreateFigure', false); % color can be adjusted by modifying the 'FaceColor', 'EdgeColor' and 'Color' parameters in "../util_functions/GU_stairsHistPlot.m"
xlim([-1, 80]);
xticks(0:10:80);
yyaxis left
yticks(0:0.01:0.05);
pbaspect([1,1,1]);
set(gca, 'Position', [.14,.12,.75,.75]);
box on
save_filename = ['Density_Outer_' save_str  '.svg'];
saveas(gcf, [resultsdir save_filename]);

%% Density Scatter Plot for each video (39 videos in total)

fig = figure();
fig.Position = [1000, 1000, 1020, 1320];
for i = 1:numel(csv_filenames)
    sp = subplot(8, 5, i);

    % Adjust subplot size
    if i <= 10
        if mod(i, 5) == 0
            sp.Position = sp.Position + [0 0.0075 0.0375 0.0375];
        else
            sp.Position = sp.Position + [0 0 0.05 0.05];
        end
    else
        if mod(i, 5) == 0
            sp.Position = sp.Position + [0 -0.0005 0.05 0.05];
        else
            sp.Position = sp.Position + [0 -0.005 0.06 0.06];
        end
    end
    densityScatterChart(double(per_exp_pair{i,1}(1,1:end)), double(per_exp_pair{i,1}(2,1:end)), 'Colormap', jet, 'ColorbarVisible', 0); % marker size can be adjusted by modifying the 'SizeData' parameter in "../util_functions/densityScatterChart.m"
    xlim([-1, 80]);
    ylim([-1, 80]);

    % Axis ticks only for the subplots along the edges
    if mod(i, 5) == 1
        yticks(0:40:80);
    else
        yticks(0:0:0);
    end

    if i > 34
        xticks(0:40:80);
    else
        xticks(0:0:0);
    end

    clim([0 5.5e-3]);
    set(gca, 'PlotBoxAspectRatio', [1,1,1]);
end
h = axes(fig,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
ylabel(h,'Distance from Outer Edge (\mum)', 'FontName', 'Arial');
xlabel(h,'Distance from Lumen (\mum)', 'FontName', 'Arial');
title(h,'Density Scatter Plot (Bapta 140 \mum)', 'FontName', 'Arial', 'Units', 'normalized', 'Position', [0.5, 1.05, 0], 'FontWeight', 'Normal');
c = colorbar(h, 'Position', [0.93 0.22 0.01 0.65], 'Ticks', 0:0.5e-3:5.5e-3);
colormap(c,'jet')
clim(h,[0 5.5e-3]);
save_filename = ['Density_Scatter_individual_1-39_' save_str  '.svg'];
saveas(gcf, [resultsdir save_filename]);