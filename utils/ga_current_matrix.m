%% GA_3D_current_matrix: computes 2D surface current trajectories using Ariane
% Create structure containing: 	.curr (3D matrix of dimensions lat,lon,time; complex numbers with real=u, imag=v)
%										.lat, .lon, .time (lon, lat, time are 1D vectors)
% Input of ga_write_ariane_current.m 


clear all;
close all;

disp('transfo current ...')


%% import data 
% https://tds.aviso.altimetry.fr/thredds/catalog/dataset-l3-swot-karin-nadir-v0_3-1d-expert/catalog.html

file_curr = 'SWOT_L3_LR_SSH_Expert_506_003_20230429T184429_20230429T193535_v0.3.nc';
%ncinfo(file_curr)
%cddisp(file_curr)


%% extracting variables

lon = double(ncread(file_curr,'longitude')) ;
lat = double(ncread(file_curr,'latitude')) ;

time_curr = double(ncread(file_curr, 'time')) ;
time_curr = time_curr+datenum(2023,01,01);

% check time: convert it into datetime
base_date = datetime('1950-01-01 00:00:00.0', 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS') ;
time_curr_cal = base_date + seconds(time_curr * 1e-6) ; % convert microseconds into seconds

%%

%to avoid repetition of lon
lon = lon(lon>=0);
lon_ok = lon(1:2:length(lon));
%lon_ok = lon_ok(lon_ok <=7);
lat_ok = lat(1:2:length(lat));
%lat_ok = lat_ok(lat_ok >= 36 & lat_ok <=40);
iter_lat=[];
for i = 1:length(lat_ok)
    if lat_ok(i) >= 36 && lat_ok(i) <=40
        iter_lat(length(iter_lat)+1)=i;
    end
end
iter_lon=[];
for i = 1:length(lon_ok)
    if lon_ok(i) >= 0 && lon_ok(i) <= 7
        iter_lon(length(iter_lon)+1)=i;
    end
end
% lat = lat_ok(lat_ok>=36.02 & lat_ok<=39.98);
% lon = lon_ok(lon_ok>=0 & lon_ok<=6.98);

u0 = ncread(file_curr,'ugos') ;
u0 = u0(iter_lon,iter_lat,120:140);
u0 = permute(u0,[2 1 3]);
v0= ncread(file_curr,'vgos') ;
v0 = v0(iter_lon,iter_lat,120:140);
v0 = permute(v0,[2 1 3]);


%%%%%% FORMATING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curr_struc=struct('lat',lat_ok,'lon',lon_ok,'time',time_curr,'curr',complex(u0,v0));
curr_struc = permute(curr_struc,[2 1]);
save('../inputs/curr_struc')
disp('transfo current : done')