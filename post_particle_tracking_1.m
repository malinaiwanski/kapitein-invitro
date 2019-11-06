
%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MKI and CC 2019-11-05

%%
% TO DO:
% filter out trajectories not along a MT

%% This code reads in output from the DoM Utrecht plug-in for detection and tracking of particles for a single movie.
% It reads in all the data and re-organizes it into a format used by cumulative_track_analysis_2.
% It can be used to visualize information about the trajectories and trial parameters.
% Filtering done in this code can be controlled using Options and Parameters

% To use this file you will need the following files:
% Results from DoM Utrecht particle localization and tracking (csv)
% MT x,y positions
%

clear all, close all
%% Options (make 0 to NOT perform related action, 1 to perform)
zplot = 1;
zsave = 0;

filt_cross_mt = 1;

set(0,'DefaultFigureWindowStyle','docked')
%addpath(' ')

%% Parameters
pixel_size = 64.0; %pixel size of camera [nm]
exp_time = 0.1; %exposure time [s]
num_frames = 600; %number of frames in movie [frames]
num_pix_x = 512;
num_pix_y = 512;

min_duration = 3; %minimum length of track to be analyzed [frames]
end_dist = 200; %maximum distance to first/last point of a MT for a spot localization to be considered to be at MT end [nm]
mt_cross_dist = 200; %distance between points on two MTs for them to be considered too close/crossing --> filtered out [nm]

%% Movie to analyze
motor = 'kif1a';
mt_type = '1cycle_cpp';
date = '2019-10-30';
filenum = 1;

%% Load data
dirname =strcat('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data','\',date,'\',motor,'\',mt_type,'\');

% Read in microtubule data
% mt_file = dir(fullfile(dirname,'MT*.txt')); %finds appropriate file
% fid=fopen(fullfile(dirname,mt_file(filenum).name)); %opens the specified file in the list and imports data
% temp_mt_data = textscan(fid,'%s %s %s %s','HeaderLines',1,'Delimiter',',','EndOfLine','\n','CommentStyle','C2'); %cell with columns %1=id %2=roi_name %3=x %4=y
% fclose(fid);
% mts = filter_mts(temp_mt_data, filt_cross_mt, mt_cross_dist, pixel_size, zplot);

% Read in particle data
motor_file = dir(fullfile(dirname,'DoM*.csv')); %finds appropriate files
fid=fopen(fullfile(dirname,motor_file(filenum).name)); %opens the specified file in the list and imports data

%for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% motor_file = dir(fullfile('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data\2019-10-30\kif1a\1cycle_cpp','DoM*.csv'));
% fid=fopen(fullfile(motor_file(1).name));

temp_motor_dat = textscan(fid,'%s %s %s %s %s %s %s %*s %*s %s %s %s %s %s %s %s %s %s %s %s %s %*s %s %s %s %*s','HeaderLines',1,'Delimiter',',_','EndOfLine','\n'); 
fclose(fid);
temp_motor_dat = reshape([temp_motor_dat{:}],size(temp_motor_dat{1,1},1),size(temp_motor_dat,2));
motor_dat = zeros(size(temp_motor_dat,1),size(temp_motor_dat,2));
for i = 1:size(temp_motor_dat,1)
    for j = 1:size(temp_motor_dat,2)
        motor_dat(i,j) = str2double(temp_motor_dat{i,j});
    end
end
% array of doubles with columns: %1=x[px] %2=y[px] %3=frame_num %4=x[nm] %5=x_loc_error[nm] %6=y[nm] %7=y_loc_error[nm] %8=amplitude_fit %9=amplitude_error %10=BG_fit %11=BG_error %12=sd_x[nm] %13=sd_x_error[nm] %14=sd_y[nm] %15=sd_y_error[nm] %16=false_positive %17=integrated_intens %18=SNR %19=R2_fit %20=track_ID %21=particle_ID %22=track_length
ntracks = motor_dat(end,20); %number of tracks
cmap=colormap(lines(ntracks));

%% Initialize figures
if zplot ~= 0
%     figure, mt_plot=gcf;
%     xlabel('x (nm)'), ylabel('y (nm)'), title('Microtubules')
%     xlim([0 num_pix_x*pixel_size]), ylim([0 num_pix_y*pixel_size])
    
    figure, traj_plot=gcf;
    xlabel('x (nm)'), ylabel('y (nm)'), title('Trajectories')
%     xlim([0 num_pix_x*pixel_size]), ylim([0 num_pix_y*pixel_size])
    
    figure, kymo_plot=gcf;
    xlabel('time (s)'), ylabel('position (nm)'), title('Kymographs')
%     xlim([0 num_frames*exp_time]), ylim([0 inf])
end
%% Initialize variables

%% Analyze
for tk = 1:ntracks
    skip_tk = 0; % 1 if track fails to meet 1 or more of the filtering conditions (incl. min duratiom, present in first/last frame)
    cens_tk = 0; %1 if spots in track too close to MT end (i.e. run length is limited by MT length -> censored)
    ind_tk = find(motor_dat(:,20)==tk); %indices of track k within motor_dat
    spots_tk = motor_dat(ind_tk,21); %particle IDs in track k
    frame_tk = motor_dat(ind_tk,3); %frames of track k
    x_tk = motor_dat(ind_tk,4); %x localizations of particles in track k [nm]
    y_tk = motor_dat(ind_tk,6); %y localizations of particles in track k [nm]
    amp_tk = motor_dat(ind_tk,8); %amplitude fit of particles in track k
    intintens_tk = motor_dat(ind_tk,17); %integrated intensity of particles in track k
    duration_tk = motor_dat(ind_tk(1),22); %duration of track k [frames]
    
    % FILTER: do not analyze track if too short
    if duration_tk < min_duration
        skip_tk = 1;
        continue
    end
    
    % FILTER: do not analyze track if present in first or last frame of movie (censored)
    if frame_tk(1) == 1 || frame_tk(end) == num_frames
        skip_tk = 1;
        continue
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FILTER: do not analyze track if not near a MT - choose MT along which the given track is, skip track if not near a MT
    % assign which index within mts as mt_tk so that mts{mt_tk} gives right MT
    % assign x,y positions and mt_coords
    
    % for testing, fit "MT" to trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pseudo_mt = polyfit(x_tk,y_tk,1);
    mt_ycoords = polyval(pseudo_mt,x_tk);
    mt_coords = [x_tk, mt_ycoords];
    [~,uni_ind] = unique(mt_coords(:,1),'stable'); %find repeated x-values in MT
    uni_ind = sort(uni_ind);
    mt_coords = mt_coords(uni_ind,:);
    
    % Find points of track considered to be at end of MT (no distinction made between plus and minus end)
    motorsq = [x_tk,y_tk];
    mtendsq = [mt_coords(1,:);mt_coords(end,:)];
    [points_mtend,dist_mtend] = rangesearch(motorsq,mtendsq,end_dist,'Distance','euclidean'); %points within specified nm of either end of MT (cell with entry for each end) %distance to MT end (cell with entry for each end)
    if size(points_mtend{1,1},2) ~= 0 || size(points_mtend{2,1},2) ~= 0 %there are points near at least one of the ends
        ind_mtend = [sort(points_mtend{1,1}), sort(points_mtend{2,1})]; %index of motorsq near MT end
        loc_mtend = [motorsq(ind_mtend,1), motorsq(ind_mtend,2)]; %x,y values of spots in motorsq near MT end
        cens_tk = 1;
    else
        ind_mtend = [];
        loc_mtend = [];
        cens_tk = 0;
    end
    all_ind = 1:1:duration_tk; %all indices (of tk variables)
    ind_alongmt = setdiff(all_ind,ind_mtend,'stable'); %indices of track (within tk variables) along length of MT

    % project track along MT
    fit_vector = [];
    motor_vector = [0,0;diff(x_tk),diff(y_tk)];
    mt_vector = [0,0;diff(mt_coords)];
    % polynomial fit to vector (gives MT coordinates for each spot localized in track)
    fit_vector(:,1) = interp1(mt_coords(:,1),mt_vector(:,1),x_tk,'pchip'); 
    fit_vector(:,2) = interp1(mt_coords(:,1),mt_vector(:,2),x_tk,'pchip');
    fit_norm = repmat(sqrt(sum(fit_vector.^2,2)),1,2);
    fit_vector_norm = fit_vector./fit_norm; %normalize
    %project
    delp = zeros(1,numel(x_tk));
    delp_off = zeros(1,numel(x_tk));
    for kv = 1:numel(x_tk)
        delp(kv) = motor_vector(kv,:)*fit_vector_norm(kv,:)'; %displacement projected onto MT axis
        delp_off(kv) = motor_vector(kv,:)*(fit_vector_norm(kv,:)*[0,1;-1,0])'; %displacement orthogonal to MT axis
    end
    position = cumsum(delp); %position of track along MT vector
    position_off = cumsum(delp_off); %position of track perpendicular to MT vector
    
    % Plot trajectory
    if zplot ~= 0
        figure(traj_plot), hold on, plot(x_tk,y_tk,'.-','Color',cmap(tk,:)) %plots trajectories
        figure(kymo_plot), hold on, plot(frame_tk.*exp_time,position,'.-','Color',cmap(tk,:)) %plots "kymographs"
    end
    
end

%% Save data
if zsave ~= 0
    
end