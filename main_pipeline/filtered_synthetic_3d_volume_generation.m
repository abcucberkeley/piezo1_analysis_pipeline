%% 3d synthetic volume generation pipeline for puncta detections

% Sayan Seal, Jan 2024

% Add necessary paths; the following paths are relative to the current directory being the "main_pipeline/" subdirectory
addpath(genpath('../LLSM5DTools/'));
addpath(genpath('../util_functions/'));

%% Load data

rt = '/clusterfs_m/nvme/sayan/20220601_Korra_MedhaGroupVisit/volumetric_2023_1/exp1_140um_2/HaloLigand/results/'; % directory containing the results obtained from point detection performed using LLSM5DTools

data = load([rt 'Detection3D.mat']).frameInfo;

%% Extract the necessary information required for filtering and creating the synthetic volumes

x = data.x;
y = data.y;
z = data.z;
A = data.A;
c = data.c;
s = data.s;
psf = data.isPSF;

%% Extract dimensions of the image stack

sz = getImageSize('/clusterfs_m/nvme/sayan/20220601_Korra_MedhaGroupVisit/volumetric_2023/exp1_140um_2/HaloLigand/DS/v_CamB_ch0_CAM1_stack0000_642nm_0000000msec_0002503655msecAbs_-01x_-01y_-01z_0001t.tif');

%% Filtering PIEZO1 puncta detections

cond_piezo1 = find(c>=100 & c<=155 & A<=65 & psf==1);
x_piezo1 = x(cond_piezo1);
y_piezo1 = y(cond_piezo1);
z_piezo1 = z(cond_piezo1);
A_piezo1 = A(cond_piezo1);
psf_piezo1 = psf(:,cond_piezo1);

% Create synthetic volume with predicted intensity profile
syn_vol_piezo1 = zeros(sz);
count = 0;
tic
for i = 1:numel(x_piezo1)
    sig_xy = s(1);
    sig_z = s(2);
    gf_sz = ceil(max(sig_xy, sig_z)*6);
    if mod(gf_sz, 10)
        gf_sz = gf_sz + (10 - mod(gf_sz, 10));
    end
    gf_sz = gf_sz + 1;
    g_vol = zeros(gf_sz,gf_sz,gf_sz); % volume to get intensity profile in a small neighborhood of the detection
    mid_pt = floor(gf_sz/2)+1;

    % Use the sigma values and the fitted amplitude of the intensity to place the Gaussian integral at the center, and then use a Gaussian filter to retrieve the
    % amplitude at the center along the intensity values in the local neighborhood
    g_vol(mid_pt,mid_pt,mid_pt) = GU_calcGaussianIntegral(A_piezo1(i), [sig_xy;sig_xy;sig_z]);
    gf_vol = filterGauss3D(g_vol, [sig_xy, sig_xy, sig_z]);

    % Ensure that the neighborhood is within the image dimensions
    if round(y_piezo1(i))-(mid_pt-1) >= 1 && round(y_piezo1(i))+(mid_pt-1) <= sz(1) && round(x_piezo1(i))-(mid_pt-1) >= 1 && round(x_piezo1(i))+(mid_pt-1) <= sz(2) && round(z_piezo1(i))-(mid_pt-1) >= 1 && round(z_piezo1(i))+(mid_pt-1) <= sz(3)
        count = count + 1;
        syn_vol_piezo1(round(y_piezo1(i))-(mid_pt-1):round(y_piezo1(i))+(mid_pt-1), round(x_piezo1(i))-(mid_pt-1):round(x_piezo1(i))+(mid_pt-1), round(z_piezo1(i))-(mid_pt-1):round(z_piezo1(i))+(mid_pt-1)) = syn_vol_piezo1(round(y_piezo1(i))-(mid_pt-1):round(y_piezo1(i))+(mid_pt-1), round(x_piezo1(i))-(mid_pt-1):round(x_piezo1(i))+(mid_pt-1), round(z_piezo1(i))-(mid_pt-1):round(z_piezo1(i))+(mid_pt-1)) + gf_vol;
    end
end
toc

% Dilate using a 3d kernel if needed
% SE = strel("sphere",3); 
% syn_vol_piezo1_dil = imdilate(syn_vol_piezo1, SE);

writetiff(uint16(syn_vol_piezo1), [rt 'synthetic_vol_piezo1_puncta_detections.tif']);

%% Extracting autofluorescence blobs

cond_blob = find(~(c>=100 & c<=155 & A<=65 & psf==1)); % not of cond_piezo1 
x_blob = x(cond_blob);
y_blob = y(cond_blob);
z_blob = z(cond_blob);
A_blob = A(cond_blob);
psf_blob = psf(:,cond_blob);

% Create synthetic volume with predicted intensity profile
syn_vol_blob = zeros(sz);
count = 0;
tic
for i = 1:numel(x_blob)
    sig_xy = s(1);
    sig_z = s(2);
    gf_sz = ceil(max(sig_xy, sig_z)*6);
    if mod(gf_sz, 10)
        gf_sz = gf_sz + (10 - mod(gf_sz, 10));
    end
    gf_sz = gf_sz + 1;
    g_vol = zeros(gf_sz,gf_sz,gf_sz); % volume to get intensity profile in a small neighborhood of the detection
    mid_pt = floor(gf_sz/2)+1;

    % Use the sigma values and the fitted amplitude of the intensity to place the Gaussian integral at the center, and then use a Gaussian filter to retrieve the
    % amplitude at the center along the intensity values in the local neighborhood
    g_vol(mid_pt,mid_pt,mid_pt) = GU_calcGaussianIntegral(A_blob(i), [sig_xy;sig_xy;sig_z]);
    gf_vol = filterGauss3D(g_vol, [sig_xy, sig_xy, sig_z]);

    % Ensure that the neighborhood is within the image dimensions
    if round(y_blob(i))-(mid_pt-1) >= 1 && round(y_blob(i))+(mid_pt-1) <= sz(1) && round(x_blob(i))-(mid_pt-1) >= 1 && round(x_blob(i))+(mid_pt-1) <= sz(2) && round(z_blob(i))-(mid_pt-1) >= 1 && round(z_blob(i))+(mid_pt-1) <= sz(3)
        count = count + 1;
        syn_vol_blob(round(y_blob(i))-(mid_pt-1):round(y_blob(i))+(mid_pt-1), round(x_blob(i))-(mid_pt-1):round(x_blob(i))+(mid_pt-1), round(z_blob(i))-(mid_pt-1):round(z_blob(i))+(mid_pt-1)) = syn_vol_blob(round(y_blob(i))-(mid_pt-1):round(y_blob(i))+(mid_pt-1), round(x_blob(i))-(mid_pt-1):round(x_blob(i))+(mid_pt-1), round(z_blob(i))-(mid_pt-1):round(z_blob(i))+(mid_pt-1)) + gf_vol;
    end
end
toc

% Dilate using a 3d kernel if needed
% SE = strel("sphere",3); 
% syn_vol_blob_dil = imdilate(syn_vol_blob, SE);

writetiff(uint16(syn_vol_blob), [rt 'synthetic_vol_autofluorescence_blobs.tif']);

%% Local density map of PIEZO1 detections

px = .108; %in um
V = 200; % in um3 search volume
r = (V*3/4/pi)^(1/3); % in um to give a V um3 spherical volume
bandwidth2 = r/px; % convert to pixel space
zAniso = .200 / .108;

filteredData_Zcorr(:,1) = y;
filteredData_Zcorr(:,2) = x;
filteredData_Zcorr(:,3) = z * zAniso;

% Use KDTree to get the absolute density in the neighboring search space
MdlKDT = KDTreeSearcher(filteredData_Zcorr);
IdxKDT = rangesearch(MdlKDT,filteredData_Zcorr,bandwidth2); % search radius = bandwidth2
filteredsynapseDensity = cellfun(@numel, IdxKDT);

syn_vol_density = zeros(sz);
tic
for i = 1:numel(x_piezo1)
    syn_vol_density(round(y_piezo1(i)), round(x_piezo1(i)), round(z_piezo1(i))) = filteredsynapseDensity(i);
end
toc

% Dilate using a 3d kernel if needed
SE = strel("sphere",3); 
syn_vol_density_dil = imdilate(syn_vol_density, SE);

writetiff(uint16(syn_vol_density_dil), [rt 'synthetic_vol_piezo1_local_density_map.tif']);

%% Density Scatter Plot of fitted amplitude and local background
% 
% fig = figure(); 
% fig.Position = [1000, 1000, 1200, 1200];
% densityScatterChart(double(A), double(c), 'Colormap', jet); % marker size
% can be adjusted by modifying the 'SizeData' parameter in "../util_functions/densityScatterChart.m"
% xlim([0, 400]);
% ylim([0, 400]);
% xlabel('Amplitude');
% ylabel('Background');
% set(gca, 'PlotBoxAspectRatio', [1,1,1]);
% set(gca, 'FontName', 'Arial', 'FontSize', 23);
% set(gca, 'Position', [.12,.15,.75,.75]);