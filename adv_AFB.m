%% Test passive process: "Corridor" hypothesis

close all, clear all

%% --- Configuration --- %%

name_curr = 'test'; % Name of current files for Ariane (from utils/transfo_current.m)

dx_init = 0.03;     % Initial meshgrid
dt = 0.1;           % Time step for Ariane
nb_days_advec = 9;  % Advection time for Ariane
dx_end = 0.1;       % Final meshgrid for maps 

t0_range = datenum(2023,4,28):datenum(2023,5,5); % Initialisation of Ariane with different t0

% Geographical limits
min_lon = 4; max_lon = 6;
min_lat = 40; max_lat = 42;

% ADT (Absolute Dynamic Topography) threshold for the front 
adt_threshold_A = -0.06;
adt_threshold_B = -0.04; 

% Date of final maps
target_date = datetime(2023,5,6);

% Variable outputs
fields = {'biomass_E_pos','biomass_E_neg','biomass_T1','biomass_T2','S', 'Region'};

% Cumulative structure with all trajectories
all_lon = [];
all_lat = [];
for f = fields
    all_val.(f{1}) = [];
end

% Set up of final maps
lon_edges = min_lon-dx_end/2:dx_end:max_lon+dx_end/2;
lat_edges = min_lat-dx_end/2:dx_end:max_lat+dx_end/2;
Longitude = min_lon:dx_end:max_lon;
Latitude = min_lat:dx_end:max_lat;
[LonGrid, LatGrid] = meshgrid(Longitude, Latitude);

% Output 2D structure
output2D = struct();
output2D.Longitude = Longitude;
output2D.Latitude = Latitude;
for f = fields
    output2D.(f{1}) = [];
end

% Ariane
global dir_ariane_global
dir_ariane_global = 'Ariane_workplace/';

%% --- Upload data --- %%

% Current for Ariane
load('inputs/curr_struc')

% Cytometry
cyto_stations = readtable('inputs/DATA_insitu/cyto_stations_AFB.csv');
E_pos = {'HfNano_biom', 'HsNano_biom', 'HflrPico_biom'};
E_neg = {'OraPicoProk_biom'};
T1 = {'RedNano_biom', 'HflrNano_biom'};
T2 = {'RedPico_biom'};
cyto_stations.E_pos = sum(cyto_stations{:, E_pos}, 2);
cyto_stations.E_neg = sum(cyto_stations{:, E_neg}, 2);
cyto_stations.T1 = sum(cyto_stations{:, T1}, 2);
cyto_stations.T2 = sum(cyto_stations{:, T2}, 2);
mean_vals = @(var, reg) mean(cyto_stations{strcmp(cyto_stations.Region, reg), var});
biomass.A2 = struct('E_pos', mean_vals('E_pos','A2'), ...
                    'E_neg', mean_vals('E_neg','A2'), ...
                    'T1', mean_vals('T1','A2'), ...
                    'T2', mean_vals('T2','A2'));
biomass.B2 = struct('E_pos', mean_vals('E_pos','B2'), ...
                    'E_neg', mean_vals('E_neg','B2'), ...
                    'T1', mean_vals('T1','B2'), ...
                    'T2', mean_vals('T2','B2'));
biomass.F2 = struct( ...
    'E_pos', mean_vals('E_pos','F2'), ...
    'E_neg', mean_vals('E_neg','F2'), ...
    'T1',    mean_vals('T1','F2'), ...
    'T2',    mean_vals('T2','F2') ...
);

% TSG salinity
tsg_transect = readtable('inputs/DATA_insitu/tsg_transect_AFB.csv');
tsg_stations = readtable('inputs/DATA_insitu/tsg_stations_AFB.csv');
S_A2 = mean(tsg_stations.S(strcmp(tsg_stations.Region, 'A2')));
S_B2 = mean(tsg_stations.S(strcmp(tsg_stations.Region, 'B2')));
S_F2 = mean(tsg_stations.S(strcmp(tsg_stations.Region, 'F2')));

%% --- Loop on t0: Advection of tracers along lagrangian trajectories --- %%

for t0 = t0_range
    
    % Find ADT file for each t0
    t0_str = datestr(t0, 'yyyymmdd'); 

    files = dir(['Ariane_workplace/currents_data/SWOT/dt_*', t0_str,'_*.nc']);
    if isempty(files), warning(['no file for ', t0_str]), continue, end
    filename = fullfile(files(1).folder, files(1).name);
    lon = ncread(filename, 'longitude');
    lat = ncread(filename, 'latitude');

    % Find the area of interest
    lon_idx = find(lon >= min_lon & lon <= max_lon);
    lat_idx = find(lat >= min_lat & lat <= max_lat);
    adt = ncread(filename, 'adt', [lon_idx(1), lat_idx(1), 1], [length(lon_idx), length(lat_idx), 1]);
    lon_sub = lon(lon_idx);
    lat_sub = lat(lat_idx);
    adt_slice = adt(:,:,1);

    % Init grid
    xvalue = min_lon:dx_init:max_lon;
    yvalue = min_lat:dx_init:max_lat;
    [X,Y] = meshgrid(xvalue, yvalue);
    point_1 = [X(:), Y(:)];
    
    % Run Ariane
    positions = ga_advection_ariane(point_1, name_curr, 'dt', dt, 'time0', t0, 'nbdays_advec', nb_days_advec);
    positions.time2D = repmat(positions.time', size(point_1,1), 1);

    N = size(point_1, 1); % Number of numerical partcile

    [biomass_E_pos, biomass_E_neg, biomass_T1, biomass_T2, S, Region] = deal(NaN(N, 1));

    for i = 1:N

        % Find the closest ADT value from the particle position
        lon_p = point_1(i, 1);
        lat_p = point_1(i, 2);
        [~, i_lon] = min(abs(lon_sub - lon_p));
        [~, i_lat] = min(abs(lat_sub - lat_p));
        val = adt_slice(i_lon, i_lat);

        if isnan(val), continue, end

        % Water mass A
        if val < adt_threshold_A
            biomass_E_pos(i) = biomass.A2.E_pos;
            biomass_E_neg(i) = biomass.A2.E_neg;
            biomass_T1(i) = biomass.A2.T1;
            biomass_T2(i) = biomass.A2.T2;
            S(i) = S_A2;
            Region(i) = 1;

        % Front F
        elseif val > adt_threshold_A && val < adt_threshold_B && lat_p < 41.3 && lon_p < 5 % F %&& lon_p < 4.5 
            biomass_E_pos(i) = biomass.F2.E_pos;
            biomass_E_neg(i) = biomass.F2.E_neg;
            biomass_T1(i) = biomass.F2.T1;
            biomass_T2(i) = biomass.F2.T2;
            S(i) = S_F2;
            Region(i) = 2;
        
        % Latitude criterion to avoid the small scale structure detected in the water mass A
        elseif val > adt_threshold_A && val < adt_threshold_B && lat_p > 41.3 
            biomass_E_pos(i) = biomass.A2.E_pos;
            biomass_E_neg(i) = biomass.A2.E_neg;
            biomass_T1(i) = biomass.A2.T1;
            biomass_T2(i) = biomass.A2.T2;
            S(i) = S_A2;
            Region(i) = 1;

        % Water mass B
        elseif val  > adt_threshold_B
            biomass_E_pos(i) = biomass.B2.E_pos;
            biomass_E_neg(i) = biomass.B2.E_neg;
            biomass_T1(i) = biomass.B2.T1;
            biomass_T2(i) = biomass.B2.T2;
            S(i) = S_B2;
            Region(i) = 3;
        end
    end

    % Add tracers to the position file
    for f = fields
        positions.(f{1}) = repmat(eval(f{1}), 1, size(positions.lon2D, 2));
    end

    % Extract the positions of the target date
    time_vec = dateshift(datetime(positions.time2D, 'ConvertFrom', 'datenum'), 'start', 'day');
    idx_6mai = time_vec == target_date;
    if ~any(idx_6mai), continue, end

    all_lon = [all_lon; positions.lon2D(idx_6mai)];
    all_lat = [all_lat; positions.lat2D(idx_6mai)];
    for f = fields
        all_val.(f{1}) = [all_val.(f{1}); positions.(f{1})(idx_6mai)];
    end
end

% Save
all_val.Longitude = all_lon;
all_val.Latitude = all_lat;
save('outputs/positions_AFB.mat', 'all_val');

%% --- Gridded map with histcn --- %%

for f = fields
    H = histcn([all_val.Longitude all_val.Latitude], lon_edges, lat_edges, ...
               'AccumData', all_val.(f{1}), 'Fun', @(x) mean(x, 'omitnan'));
    H(H == 0) = NaN;
    output2D.(f{1}) = H; % 2D matrix (lon*lat)
end

%%  Add Region labels

% Count the number of A and B
regionA = histcn([all_val.Longitude all_val.Latitude], lon_edges, lat_edges, 'AccumData', all_val.Region==1);
regionF = histcn([all_val.Longitude all_val.Latitude], lon_edges, lat_edges, 'AccumData', all_val.Region==2);
regionB = histcn([all_val.Longitude all_val.Latitude], lon_edges, lat_edges, 'AccumData', all_val.Region==3);

% Proportions
denom = regionA+regionF+regionB;
output2D.percentA = regionA ./ denom;
output2D.percentF = regionF ./ denom;
output2D.percentB = regionB ./ denom;

% Add flag for A, F, B (Find pixels where A, F or B is the most represented label
percent_stack = cat(3, output2D.percentA, output2D.percentF, output2D.percentB);
percent_stack(isnan(percent_stack)) = -Inf;
[val_max, idx_max] = max(percent_stack, [], 3);
idx_max(val_max == -Inf) = NaN;
Region_flag = nan(size(output2D.percentA));
Region_flag(idx_max == 1) = 1;
Region_flag(idx_max == 2) = 2;
Region_flag(idx_max == 3) = 3;
output2D.Region_flag = Region_flag;

save('outputs/output2D_AFB.mat', 'output2D');

%% Figures

figure;
hold on;

% for i = 1:N
%     scatter(positions.lon2D(i, :), positions.lat2D(i, :), 0.5, positions.Region(i, :), 'filled');
% end

% Add ship track
colors = zeros(length(tsg_transect.Region), 3); 
colors(strcmp(tsg_transect.Region, 'A'), :) = repmat([0.3010 0.7450 0.9330], sum(strcmp(tsg_transect.Region, 'A')), 1);
colors(strcmp(tsg_transect.Region, 'F'), :) = repmat([0.8500 0.3250 0.0980], sum(strcmp(tsg_transect.Region, 'F')), 1);
colors(strcmp(tsg_transect.Region, 'B'), :) = repmat([0.4660 0.6740 0.1880], sum(strcmp(tsg_transect.Region, 'B')), 1);
scatter(tsg_transect.Longitude, tsg_transect.Latitude, 5, colors, 'filled');
colors_stations = zeros(length(tsg_stations.Region), 3);
colors_stations(strcmp(tsg_stations.Region, 'A2'), :) = repmat([0.3010 0.7450 0.9330], sum(strcmp(tsg_stations.Region, 'A2')), 1);
colors_stations(strcmp(tsg_stations.Region, 'F2'), :) = repmat([0.8500 0.3250 0.0980], sum(strcmp(tsg_stations.Region, 'F2')), 1);
colors_stations(strcmp(tsg_stations.Region, 'B2'), :) = repmat([0.4660 0.6740 0.1880], sum(strcmp(tsg_stations.Region, 'B2')), 1);
scatter(tsg_stations.Longitude, tsg_stations.Latitude, 5, colors_stations, 'filled');

cmap1 = colormap([
    0.6 0.8 1.0
    1.0 0.6 0.6
    0.6 0.9 0.6
]);
colorbar('Ticks', [1 2 3], 'TickLabels', {'A', 'F', 'B'});

title('Trajectories');

print('-djpeg','-r300','figures/adv_AFB_traj.jpg')

fields = {'biomass_E_pos','biomass_E_neg','biomass_T1','biomass_T2','S'};

for f = fields
    figure;
    pcolor(LonGrid, LatGrid, output2D.(f{1})'); shading flat;
    colorbar; title(['6 mai - ', f{1}, ' (mean)']);
    print('-djpeg','-r300',['figures/composite_mean_', f{1}, '.jpg']);
end

% Carte des labels
figure;
Region_num = double(output2D.Region_flag)';
pcolor(LonGrid, LatGrid, Region_num);
shading flat;
colormap(cmap1);
colorbar('Ticks', [1 2 3], 'TickLabels', {'A','F','B'});
title('6 mai - Label');
xlabel('Longitude');
ylabel('Latitude');
hold on;
scatter(tsg_transect.Longitude, tsg_transect.Latitude, 5, colors, 'filled');
scatter(tsg_stations.Longitude, tsg_stations.Latitude, 5, colors_stations, 'filled');
print('-djpeg','-r300','figures/adv_label_AFB.jpg');

% Carte label end
figure;
scatter(all_val.Longitude, all_val.Latitude, 30, all_val.Region, 'filled');
xlim([min_lon max_lon]); ylim([min_lat max_lat]);
colormap(cmap1);
colorbar('Ticks', [1 2 3], 'TickLabels', {'A','F','B'});
colorbar;
title('6 mai - end position - labels');
print('-djpeg','-r300','figures/adv_AFB_finalscatter_label.jpg');