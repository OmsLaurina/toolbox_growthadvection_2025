%% transfor_currents: create current matrice for Ariane2D
%
% a 3D matrix is needed with longitude, latitude and time
% Adapted from N.Kientz
%

clear all;
close all;

disp('transfo current ...')

% Set directory where Ariane configuration files are (used in most functions). This is different from the directory where Ariane is installed.
% There must be a directory "currents_data" inside dir_ariane_global where currents netcdf files are saved (see ga_write_ariane_currents)
global dir_ariane_global
dir_ariane_global='Ariane_workplace/';

% Liste des fichiers (à adapter selon vos fichiers et votre répertoire)
file_list = dir('Ariane_workplace/currents_data/SWOT/dt*.nc');
num_files = length(file_list);

% Initialisation des variables pour stocker l'ensemble des données
lat_combined = [];
lon_combined = [];
u_all = [];
v_all = [];
time_combined = [];
time_str_combined = {};

for i = 1:num_files
    % Charger le fichier courant
    file_curr = fullfile(file_list(i).folder, file_list(i).name);
    
    % Extraction des coordonnées géographiques et de la date
    lat = double(ncread(file_curr, 'latitude'));
    lon = double(ncread(file_curr, 'longitude'));
    time_curr = double(ncread(file_curr, 'time')) + datenum(1950,01,01);
    time_curr_str = datestr(time_curr, 'dd/mm/yyyy');
    
    % Définir les indices pour les limites de la grille
    lat_i = find(lat >= 39 & lat <= 43.2);
    lon_i = find(lon >= 3 & lon <= 6);
    lat = lat(lat_i);
    lon = lon(lon_i);
    
    % Extraction des variables u et v pour cette date
    u0 = ncread(file_curr, 'ugos', [lon_i(1), lat_i(1), 1], [length(lon_i), length(lat_i), 1]);
    v0 = ncread(file_curr, 'vgos', [lon_i(1), lat_i(1), 1], [length(lon_i), length(lat_i), 1]);
    u0 = permute(u0, [2, 1]);
    v0 = permute(v0, [2, 1]);


    % Ajouter les données extraites aux matrices combinées
    if isempty(lat_combined) && isempty(lon_combined)
        % Initialisation des latitudes et longitudes
        lat_combined = lat;
        lon_combined = lon;
    end
    
    u_all = cat(3, u_all, u0);
    v_all = cat(3, v_all, v0);
    time_combined = [time_combined; time_curr];
    time_str_combined = [time_str_combined; {time_curr_str}];
    i
end

% Création de la structure finale curr_stru
curr_struc = struct('lat', lat_combined, 'lon', lon_combined, ...
                   'time', time_combined, 'time_str', {time_str_combined}, ...
                   'curr', complex(u_all, v_all));

% Sauvegarde de la structure
save('outputs/curr_struc', 'curr_struc');

%% write currents
ga_write_ariane_currents(curr_struc,'test')

disp('transfo current : done')