%% Test passive process: "Boundary" hypothesis

close all, clear all

%% --- Configuration --- %%

name_curr = 'test'; % Current files for Ariane (from utils/transfo_current.m)

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
plankton_model_outputs = {'P1','P2','Z','PO4','u1','u2','g1','g2'};

% Cumulative structure with all trajectories
all_lon = [];
all_lat = [];
all_val_NPZ = struct();
for v = 1:length(plankton_model_outputs)
    all_val_NPZ.(plankton_model_outputs{v}) = [];
end

% Set up of final maps
lon_edges = min_lon-dx_end/2:dx_end:max_lon+dx_end/2;
lat_edges = min_lat-dx_end/2:dx_end:max_lat+dx_end/2;
Longitude = min_lon:dx_end:max_lon;
Latitude = min_lat:dx_end:max_lat;
[LonGrid, LatGrid] = meshgrid(Longitude, Latitude);

% Init grid for NPZ output
plankton_model_outputs = {'P1','P2','Z','PO4','u1','u2','g1','g2'};
output2D_NPZ = struct();
output2D_NPZ.Longitude = Longitude;
output2D_NPZ.Latitude = Latitude;
for v = 1:length(plankton_model_outputs)
    output2D_NPZ.(plankton_model_outputs{v}) = zeros(length(Longitude), length(Latitude));
end

% Nsupply & Gmax values
Nsupply_A = 0.052; %0.048;
Nsupply_F = 0.05; %0.052;
Nsupply_B = 0.048; %0.05;

Gmax1_A = 3.89;
Gmax1_F = 3.89+1;
Gmax1_B = 3.89;
Gmax2_A = 0.43;
Gmax2_F = 0.43+1;
Gmax2_B = 0.43; 

% Import Gmax matrices (from F_GMM.py)
Gmax1 = readmatrix('outputs/Gmax1_AFB.csv');
Gmax2 = readmatrix('outputs/Gmax2_AFB.csv');

% Ariane
global dir_ariane_global
dir_ariane_global = 'Ariane_workplace/';

%% --- Upload data --- %%

% Current for Ariane
load('inputs/curr_struc')

% TSG salinity (for figures)
tsg_transect = readtable('inputs/DATA_insitu/tsg_transect_AFB.csv');
tsg_stations = readtable('inputs/DATA_insitu/tsg_stations_AFB.csv');

%% --- Loop on t0: Advection of NPZ output along lagrangian trajectories --- %%

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

    %% NPZ model

    % Create matrix for Nsupply and Gmax
    nlon = length(lon_idx);
    nlat = length(lat_idx);
    ntime = length(positions.time);

    folder = 'Ariane_workplace/currents_data/SWOT/';
    file_prefix = 'dt_global_allsat_phy_l4_';
    file_suffix = '_20240501.nc';

    Nsupply = NaN(nlon, nlat, ntime);

    for day = 0:nb_days_advec-1
        current_date = t0 + days(day);
        file_date_str = datestr(current_date, 'yyyymmdd');
        filename = fullfile(folder, [file_prefix, file_date_str, file_suffix]);
        if ~isfile(filename), continue, end

        adt_day = ncread(filename, 'adt', [lon_idx(1), lat_idx(1), 1], [nlon, nlat, 1]);

        lat_Nsupply = ncread(filename, 'latitude');
        lat_idx_Nsupply = find(lat_Nsupply >= min_lat & lat_Nsupply <= max_lat);
        lat_sub_Nsupply = lat_Nsupply(lat_idx_Nsupply);
        
        for substep = 1:(1/dt)
            t_idx = round(day/dt) + substep;
            lat_mask = repmat(lat_sub_Nsupply' > 41.5, [nlon,1]);

            Nsupply_tmp = Nsupply_A * (adt_day < adt_threshold_A) + ...
                      Nsupply_B * (adt_day > adt_threshold_B) + ...
                      Nsupply_F * (adt_day >= adt_threshold_A & adt_day <= adt_threshold_B);
            Nsupply(:, :, t_idx) = Nsupply_A .* lat_mask + Nsupply_tmp .* (~lat_mask);
        end
    end

    % Run the plankton model to get equilibrium values
    tvec = 1:dt:2000;
    output_A = ga_model_2P1Z_fromNsupplyGmax(Nsupply_A, Gmax1_A, Gmax2_A,'time', tvec);
    CI_A = [output_A.PO4(end), output_A.P1(end), output_A.P2(end), output_A.Z(end)];
  
    output_F = ga_model_2P1Z_fromNsupplyGmax(Nsupply_F, Gmax1_F, Gmax2_F, 'time', tvec);
    CI_F = [output_F.PO4(end), output_F.P1(end), output_F.P2(end), output_F.Z(end)];

    output_B = ga_model_2P1Z_fromNsupplyGmax(Nsupply_B, Gmax1_A, Gmax2_A,'time', tvec);
    CI_B = [output_B.PO4(end), output_B.P1(end), output_B.P2(end), output_B.Z(end)];

    % Init NPZ
    Npart = size(positions.lon2D, 1);
    Ntime = length(positions.time);
    results = struct();
    for v = 1:length(plankton_model_outputs)
        results.(plankton_model_outputs{v}) = NaN(Npart, Ntime);
    end

    %% Ariane + NPZ
    for i = 1:Npart

        Nsupply_series = NaN(1, Ntime);

        % Generate the time series for Nsupply and Gmax
        for t = 1:Ntime
            lon_i = positions.lon2D(i,t);
            lat_i = positions.lat2D(i,t);
            [~, idx_lon] = min(abs(lon_sub - lon_i));
            [~, idx_lat] = min(abs(lat_sub - lat_i));

            Nsupply_series(t) = Nsupply(idx_lon, idx_lat, t);
        end

        Gmax1_series = Gmax1(idx_lon, idx_lat);
        Gmax2_series = Gmax2(idx_lon, idx_lat);

        % Get the initial equilibrium concentrations for A and B
        Nsupply_ini = Nsupply_series(1);
        if Nsupply_ini == Nsupply_A
            PO4_ini = CI_A(1); P1_ini = CI_A(2); P2_ini = CI_A(3); Z_ini = CI_A(4);

        elseif Nsupply_ini == Nsupply_F
            PO4_ini = CI_F(1); P1_ini = CI_F(2); P2_ini = CI_F(3); Z_ini = CI_F(4);

        elseif Nsupply_ini == Nsupply_B
            PO4_ini = CI_B(1); P1_ini = CI_B(2); P2_ini = CI_B(3); Z_ini = CI_B(4);
        else
            continue
        end

        % Run the plankton model
        output = ga_model_2P1Z_fromNsupplyGmax(Nsupply_series, Gmax1_series, Gmax2_series,...
            'P1_ini', P1_ini, 'P2_ini', P2_ini, 'PO4_ini', PO4_ini, 'Z_ini', Z_ini, ...
            'time', positions.time - positions.time(1));

        % Stock results
        for v = 1:length(plankton_model_outputs)
            varname = plankton_model_outputs{v};
            if isfield(output, varname)
                results.(varname)(i, :) = output.(varname)(1:Ntime);
            end
        end
    end

    results.Longitude = positions.lon2D;
    results.Latitude = positions.lat2D;
    results.Time = positions.time;
    
    % Extract results of the target date
    time_vec = dateshift(datetime(positions.time2D, 'ConvertFrom', 'datenum'), 'start', 'day');
    idx_6mai = time_vec == target_date;
    if ~any(idx_6mai), continue, end
    
    all_lon = [all_lon; results.Longitude(idx_6mai)];
    all_lat = [all_lat; results.Latitude(idx_6mai)];
    for v = 1:length(plankton_model_outputs)
        varname = plankton_model_outputs{v};
        val = results.(varname)(idx_6mai);
        all_val_NPZ.(varname) = [all_val_NPZ.(varname); val];
    end
end

% Save 
all_val_NPZ.Longitude = all_lon;
all_val_NPZ.Latitude = all_lat;
save('outputs/positionsNPZ_AFB.mat', 'all_val_NPZ');

%% --- Gridded map --- %%

for v = 1:length(plankton_model_outputs)
    varname = plankton_model_outputs{v};
    H = histcn([all_lon all_lat], lon_edges, lat_edges, ...
               'AccumData', all_val_NPZ.(varname), 'Fun', @(x) mean(x, 'omitnan'));
    H(H == 0) = NaN;
    output2D_NPZ.(varname) = H;  % 2D matrix (lon*lat)
end

% Compute R-Ratio
output2D_NPZ.R_P1 = output2D_NPZ.P1 ./ (output2D_NPZ.P1 + output2D_NPZ.P2);
output2D_NPZ.R_P2 = output2D_NPZ.P2 ./ (output2D_NPZ.P1 + output2D_NPZ.P2);
output2D_NPZ.Btot = output2D_NPZ.P1 + output2D_NPZ.P2;

% Save output
save('outputs/output2D_AFB_NPZ_activereactive.mat', 'output2D_NPZ')

%% Figures 

% for ship track
colors = zeros(length(tsg_transect.Region), 3); 
colors(strcmp(tsg_transect.Region, 'A'), :) = repmat([0.3010 0.7450 0.9330], sum(strcmp(tsg_transect.Region, 'A')), 1);
colors(strcmp(tsg_transect.Region, 'F'), :) = repmat([0.8500 0.3250 0.0980], sum(strcmp(tsg_transect.Region, 'F')), 1);
colors(strcmp(tsg_transect.Region, 'B'), :) = repmat([0.4660 0.6740 0.1880], sum(strcmp(tsg_transect.Region, 'B')), 1);
colors_stations = zeros(length(tsg_stations.Region), 3);
colors_stations(strcmp(tsg_stations.Region, 'A2'), :) = repmat([0.3010 0.7450 0.9330], sum(strcmp(tsg_stations.Region, 'A2')), 1);
colors_stations(strcmp(tsg_stations.Region, 'F2'), :) = repmat([0.8500 0.3250 0.0980], sum(strcmp(tsg_stations.Region, 'F2')), 1);
colors_stations(strcmp(tsg_stations.Region, 'B2'), :) = repmat([0.4660 0.6740 0.1880], sum(strcmp(tsg_stations.Region, 'B2')), 1);

% Maps
varnames_fig = {'PO4', 'P1', 'P2', 'Z', 'R_P1', 'R_P2', 'Btot'};
titles = {'Masse PO4','Biomasse P1','Biomasse P2', 'Biomasse Z', 'R-ratio P1', 'R-ratio P2', 'Biomasse Tot P'};
filenames = {'advgrowth_biom_Nsupply_AB_PO4.jpg','advgrowth_biom_Nsupply_AB_P1.jpg',...
    'advgrowth_biom_Nsupply_AB_P2.jpg', 'advgrowth_biom_Nsupply_AB_Z.jpg',...
    'advgrowth_biom_Nsupply_AB_RP1.jpg', 'advgrowth_biom_Nsupply_AB_RP2.jpg', 'advgrowth_biom_Nsupply_Btot.jpg'};

for i = 1:length(varnames_fig)
    figure;
    pcolor(LonGrid, LatGrid, output2D_NPZ.(varnames_fig{i})');
    shading flat;
    colorbar;
    title(['6 mai - ', titles{i}]);
    hold on;
    scatter(tsg_transect.Longitude, tsg_transect.Latitude, 5, colors, 'filled');
    scatter(tsg_stations.Longitude, tsg_stations.Latitude, 5, colors_stations, 'filled');
    if strcmp(varnames_fig{i}, 'S')
        caxis([38.1 38.5]);
    end
    print('-djpeg','-r300',['figures/', filenames{i}]);
end