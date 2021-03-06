
%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski 2019-11-05
% Projection along MT taken from AG Hendricks Lab (McGill)

%%
% TO DO:
% in filter_mts, filter out MTs that run out of edge of FOV
% landing rate per MT with time (also for whole MT for whole movie div by L MT and L mov)

%% This code reads in output from the DoM Utrecht plug-in for detection and tracking of particles for a single movie.
% It reads in all the data and re-organizes it into a format used by cumulative_track_analysis_2.
% It can be used to visualize information about the trajectories and trial parameters.
% Filtering done in this code can be controlled using Options and Parameters

% To use this file you will need the following files:
% Results from DoM Utrecht particle localization and tracking (csv; file name must start with DoM)
% MT x,y positions (csv; file name must start with MT)
% boundaries between different sections of MT - only if zcap == 1 (csv; file name must start with int)

% You will also need to download the following functions:
% arclength from https://nl.mathworks.com/matlabcentral/fileexchange/34871-arclength
% interparc from https://nl.mathworks.com/matlabcentral/fileexchange/34874-interparc

clear all, close all
%addpath('C:\Users\6182658\OneDrive - Universiteit Utrecht\MATLAB\GitHub Codes\in-vitro-codes\kapitein-invitro') %windows
addpath('/Users/malinaiwanski/Documents/MATLAB/GitHub/kapitein-invitro') %mac
%addpath('C:\Users\6182658\OneDrive - Universiteit Utrecht\MATLAB') %windows
addpath('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/MATLAB') %mac
addpath('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/in_vitro_data') %mac
set(0,'DefaultFigureWindowStyle','docked')

%% Movies to analyze
motors = {'kif1a','kif5b'}; %
mt_types = {'cap','taxol_cap'}; %{'1cycle_cpp','2cycle_cpp','gdp_taxol'}; %
dates =  {'2019-12-09','2019-12-13','2018-10-17'}; %{'2020-02-05','2019-10-30'}; %

%%
for master_date_ind = 1:size(dates,2)
    for master_motor_ind = 1:size(motors,2)
        for master_mt_ind = 1:size(mt_types,2)

            %% Movie to analyze
            motor = motors{master_motor_ind}; %'kif5b'; %'kif1a'; %
            mt_type = mt_types{master_mt_ind}; %'gdp_taxol'; %'cap'; %'1cycle_cpp'; %'2cycle_cpp'; %
            date = dates{master_date_ind}; %'2019-10-30'; %'2019-12-09'; %

            %filenum = %5;

            %% Load data
            %dirname =strcat('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data','\',date,'\',motor,'\',mt_type,'\'); %windows
            dirname =strcat('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/in_vitro_data','/',date,'/',motor,'/',mt_type,'/'); %mac


            motor_files = dir(fullfile(dirname,'DoM*.csv')); %finds appropriate files
            num_files = size(motor_files,1);

            for filenum = 1:num_files
                clearvars -except motors mt_types dates master_date_ind master_motor_ind master_mt_ind motor mt_type date dirname filenum num_files
                close all
                
                %% Options (make 0 to NOT perform related action, 1 to perform)
                zplot = 0; %set to 1 to visualize trajectories, kymographs, etc.
                zsave = 1; %set to 1 to save the output from this file, must be done if planning to use cumulative_track_analysis_2
                zcap = 0; %set to 1 if using capped MTs
                zcurve = 0; %set to 1 if Curve Tracing was used to trace the MTs in Fiji
%                 if date == '2020-02-05'
%                     zcurve = 1;
%                 end

                % Filtering
                filt_cross_mt = 1; %ignore any tracks on MTs that are too close to another MT - set this distance in the parameters > for analysis section
                filt_short_mt = 1; %ignore any tracks on MTs that are too short - set this length in the parameters > for analysis section
                filt_rl = 1; %ignore track with a very short displacement (likely static) - set this distance in the parameters > for analysis section
                filt_tk_mt_end = 0; %ignore any tracks that are too close to the end of the MT - set this distance in the parameters > for analysis section (end_dist)

                %% Parameters
                % From imaging:
                pixel_size = 64.0; %pixel size of camera [nm]
                exp_time = 0.1; %exposure time [s]
                num_frames = 600; %number of frames in movie [frames]
                num_pix_x = 512;
                num_pix_y = 512;

                % For analysis:
                min_duration = 5; %minimum length of track to be analyzed [frames]
                min_rl = 150; %minimum run length of track to be analyzed [nm]
                end_dist = 200; %maximum distance to first/last point of a MT for a spot localization to be considered to be at MT end [nm]; if distance is < this, points in track are marked and track may be ignored if corresponding filtering parameter is set to 1
                mt_dist = 400; %maximum distance to some MT for track to considered on a MT and thus analyzed [nm]
                analyze_mt_num = -1; %specify which MT to analyze, based on MT id in text file (these start at 0); if this is -1, all MTs will be analyzed
                mt_cross_dist = 200; %distance between points on two MTs for them to be considered too close/crossing --> filtered out [nm]
                mt_min_length = 2000; %shortest MT to be analyzed, if shorter --> filtered out [nm]
                segment_boundary_dist = 200; %distance to segment boundary (i.e. between CPP seed - GDP lattice - CPP cap) to decide if/where track crosses such a boundary [nm]
                l_window = 7; %number of frames to average for sliding MSD analysis
                msd_thresh = 1.1; %alpha-value above which is processive, below which is paused
                msd_step = 0.3; %minimum threshold for findchangepts function; minimum improvement in residual error; changepoints currently based on mean values, but can use rms, std, mean and slope
                l_min = 3; %minimum distance between two changepoints - smallest duration of pause/run [frames]

                % For plotting:
                rl_binwidth = 100; %bin width for run length histograms

                %%
                % Read in microtubule data
                if zcurve == 0
                    mt_file = dir(fullfile(dirname,'MT*.csv')); %finds appropriate file
                    fid=fopen(fullfile(dirname,mt_file(filenum).name)); %opens the specified file in the list and imports data
                    temp_mt_data = textscan(fid,'%s %s %s %s','HeaderLines',1,'Delimiter',',','EndOfLine','\n','CommentStyle','C2'); %cell with columns %1=id %2=roi_name %3=x %4=y
                    fclose(fid); 
                else
                    mt_file = dir(fullfile(dirname,'MT*.csv')); %finds appropriate file
                    fid=fopen(fullfile(dirname,mt_file(filenum).name)); %opens the specified file in the list and imports data
                    temp_mt_dat = textscan(fid,'%s %*s %s %*s %s %s %*s %s %*s %*s %*s %*s','HeaderLines',1,'Delimiter',',','EndOfLine','\n');
                    fclose(fid);
                    temp_mt_dat = reshape([temp_mt_dat{:}],size(temp_mt_dat{1,1},1),size(temp_mt_dat,2));
                    response_intens{1,1} = strings(size(temp_mt_dat,1),1);
                    for i = 1: size(temp_mt_dat,1)
                        response_intens{1,1}(i) = convertCharsToStrings(temp_mt_dat{i,5}); %str2double(temp_mt_dat{i,5});
                    end
                    temp_mt_dat(:,5) = [];
                    temp_mt_dat_old = temp_mt_dat;
                    temp_mt_dat(:,1) = temp_mt_dat_old(:,2);
                    temp_mt_dat(:,2) = temp_mt_dat_old(:,1);
                    % temp_mt_data = zeros(size(temp_mt_dat,1),size(temp_mt_dat,2));
                    % for i = 1:size(temp_mt_data,1)
                    %     for j = 1:size(temp_mt_data,2)
                    %         temp_mt_data(i,j) = str2double(temp_mt_dat{i,j});
                    %     end
                    % end
                    % temp_mt_data = num2cell(temp_mt_data);
                    temp_mt_data = cell(1,4);
                    for i = 1:size(temp_mt_dat,1)
                        for j = 1:size(temp_mt_dat,2) %[1,3,4]
                            temp_mt_data{1,j} = [temp_mt_data{1,j};convertCharsToStrings(temp_mt_dat{i,j})];%str2double(temp_mt_dat{i,j})];
                        end
                    end
                end
                
                if date == '2020-02-05'
                    temp_ids = str2double(temp_mt_data{1,1}(:));
                    temp_ids = temp_ids-1;
                    temp_mt_data{1,1} = [];
                    for mti = 1:length(temp_ids)
                        temp_mt_data{1,1} = [temp_mt_data{1,1};convertCharsToStrings(num2str(temp_ids(mti)))];
                    end
                end
                
                if zcurve == 1 && zcap == 1
                    [mts, interp_mts, response_mts, skip_mts] = filter_mts([temp_mt_data,response_intens], analyze_mt_num, filt_cross_mt, mt_cross_dist, filt_short_mt, mt_min_length, pixel_size, num_pix_x, num_pix_y, zplot, zcap, zcurve);
                else
                    [mts, interp_mts, ~, skip_mts] = filter_mts(temp_mt_data, analyze_mt_num, filt_cross_mt, mt_cross_dist, filt_short_mt, mt_min_length, pixel_size, num_pix_x, num_pix_y, zplot, zcap, zcurve);
                end
                all_interp_mts = [];
                for j = 1: size(interp_mts,1)
                    all_interp_mts = [all_interp_mts;interp_mts{j}];
                end
                filt_mt_ind = find(skip_mts); % This will be changed later for capped MTs if any are found that do not have 3 or 5 segments
                num_mts = size(mts,1);
                
                % For capped MTs, read in segment data
                if zcap == 1
                    cap_file = dir(fullfile(dirname,'int*.csv')); %finds appropriate file
                    fid=fopen(fullfile(dirname,cap_file(filenum).name)); %opens the specified file in the list and imports data
                    temp_boundary_data = textscan(fid,'%s %s %s %s','HeaderLines',1,'Delimiter',',','EndOfLine','\n','CommentStyle','C2'); %cell with columns %1=id %2=roi_name %3=x %4=y
                    fclose(fid); 

                    %format xy positions
                    temp_caps_x = str2double(temp_boundary_data{1,3}(:));
                    temp_caps_y = str2double(temp_boundary_data{1,4}(:)); 
                    if zcap == 1 %most MTs vertical, flip x,y
                        temp_caps_x = str2double(temp_boundary_data{1,4}(:));
                        temp_caps_y = str2double(temp_boundary_data{1,3}(:)); 
                    end

                    temp_cap_data = [temp_caps_x,temp_caps_y];
%                     uni_caps = unique(temp_cap_data,'rows','stable'); %The ROI saver will repeat all preceding points clicked, so unique saves only the first instance of each point
                    tol = 0.001;
                    uni_caps = uniquetol(temp_cap_data,tol,'ByRows',true);
                    
                    caps_x = uni_caps(:,1).*pixel_size+1.5*pixel_size;
                    caps_y = uni_caps(:,2).*pixel_size+1.5*pixel_size;
                    cap_loc = [caps_x, caps_y];

                    %identify on which MT the identified segment boundary is
                    caps_id = [];
                    segment_boundaries = cell(num_mts,1);
                    for i = 1:length(caps_x)
                        closest_point_on_mt = zeros(num_mts,1);
                        for j = 1:num_mts
                            rep_loc = repmat(cap_loc(i,:),size(interp_mts{j},1),1);
                            diff_loc = interp_mts{j}-rep_loc;
                            cap_mt_dist = vecnorm(diff_loc,2,2);
                            closest_point_on_mt(j) = min(cap_mt_dist);
                        end
                        [minval,closemt] = min(closest_point_on_mt);
                        cap_id = closemt;
                        caps_id = [caps_id, cap_id];
                        segment_boundaries{cap_id} = [segment_boundaries{cap_id};caps_x(i),caps_y(i)];
                    end
                    
                    for i = cap_id
                        if size(segment_boundaries{i},1) ~=2 && size(segment_boundaries{i},1) ~= 4
                            disp(['MT number ', num2str(cap_id), ' does not have 3 or 5 segments. Retrace segment boundaries']);
                        skips_mts(i) = 1; % also filter out MTs if they do not have an appropriate number of segments
                        end
                    end
                    
                    filt_mt_ind = find(skip_mts);
                    
                    num_caps = length(unique(caps_id));
                    if num_caps ~= num_mts
                        disp('ERROR: Data mismatch - retrace MTs and/or MT segment boundaries')
                    end

                    boundaries_on_mt = cell(num_mts,1);
                    for j = 1:num_mts

                        %find closest point along MT to each segment boundary
                        for i = 1:size(segment_boundaries{j},1)
                            [nearpoints_mt,dist_nearpoints] = rangesearch(interp_mts{j},segment_boundaries{j},5000,'Distance','euclidean');%,'SortIndices',1); %find nearby points on MT
                            closest_mtpoint_ind = nearpoints_mt{i,1}(1);
                            closest_mtpoint = interp_mts{j}(closest_mtpoint_ind,:);
                            boundaries_on_mt{j} = [boundaries_on_mt{j}; closest_mtpoint]; 
                        end
                    end
                end
                
                %identify MT orientation and which indices in MT correspond to the MT seed
                if zcurve == 1 && zcap == 1
                    flip_this_mt = {};
                    mt_seed = {};
                    for j = 1:num_mts
                        [flip_this_mt{j}, mt_seed{j}] = find_mt_seed(interp_mts{j},boundaries_on_mt{j},response_mts{j});
                    end
                end
                
                % Read in motor concentration
                conc_file = dir(fullfile(dirname,'concentration*.csv')); %finds appropriate file
                fid=fopen(fullfile(dirname,conc_file(1).name)); %opens the specified file in the list and imports data
                temp_conc_data = textscan(fid,'%s %s','HeaderLines',1,'Delimiter',';','EndOfLine','\n','CommentStyle','C2'); %cell with columns %1=id %2=roi_name %3=x %4=y
                fclose(fid);
                conc_data = zeros(num_files,size(temp_conc_data,2));
                for i = 1:num_files %size(temp_conc_data,1)
                    for j = 1:size(temp_conc_data,2)
                        conc_data(i,j) = str2double(temp_conc_data{1,j}(i,1));
                    end
                end
                motor_concentration = conc_data(filenum,2); % [nM]
                
                % Read in particle data
                motor_file = dir(fullfile(dirname,'DoM*.csv')); %finds appropriate files
                fid=fopen(fullfile(dirname,motor_file(filenum).name)); %opens the specified file in the list and imports data

                disp('------------Reading in data from file:------------')
                disp(strcat(date,'\',motor,'\',mt_type,'\',motor_file(filenum).name))

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

                if zcap == 1
                    old_motor_dat = motor_dat;
                    motor_dat(:,1) = old_motor_dat(:,2);
                    motor_dat(:,2) = old_motor_dat(:,1);
                    motor_dat(:,4) = old_motor_dat(:,6);
                    motor_dat(:,6) = old_motor_dat(:,4);
                    motor_dat(:,5) = old_motor_dat(:,7);
                    motor_dat(:,7) = old_motor_dat(:,5);
                    motor_dat(:,12) = old_motor_dat(:,14);
                    motor_dat(:,14) = old_motor_dat(:,12);
                    motor_dat(:,13) = old_motor_dat(:,15);
                    motor_dat(:,15) = old_motor_dat(:,13);
                end

                % array of doubles with columns: %1=x[px] %2=y[px] %3=frame_num %4=x[nm] %5=x_loc_error[nm] %6=y[nm] %7=y_loc_error[nm] %8=amplitude_fit %9=amplitude_error %10=BG_fit %11=BG_error %12=sd_x[nm] %13=sd_x_error[nm] %14=sd_y[nm] %15=sd_y_error[nm] %16=false_positive %17=integrated_intens %18=SNR %19=R2_fit %20=track_ID %21=particle_ID %22=track_length
                ntracks = motor_dat(end,20); %number of tracks
                % cmap = colormap(parula(num_mts));
                % cmap = colormap(lines(num_mts));
                cmap = colormap(colorcube(num_mts));

                %% Initialize figures
                if zplot ~= 0
                    figure, traj_plot=gcf;
                    xlabel('x (nm)'), ylabel('y (nm)'), title('Trajectories')

                    figure, kymo_plot=gcf;
                    xlabel('time (s)'), ylabel('position (nm)'), title('Kymographs')

                    figure, offaxis_plot=gcf;
                    xlabel('time (s)'), ylabel('off-axis position (nm)'), title('Off-axis position')

                end
                %% Initialize variables
                traj = struct([]);
                cum_run_length = [];
                censored = []; %track reaches end of MT
                cum_censored = [];
                cum_mean_vel = [];
                cum_proc_vel = [];
                cum_pause_vel = [];
                cum_inst_vel = [];
                cum_loc_alpha = [];
                cum_proc_alpha = [];
                cum_pause_alpha = [];
                cum_association_time = [];
                cum_mts = [];
                cum_norm_landing_pos = [];
                cum_landing_dist_to_mt_end = [];
                flip_mt = zeros(1,num_mts);
                no_flip = zeros(1,num_mts);
                mt_lengths = cell(1,num_mts);
                ftk = 0; %number of tracks that passes all filtering
                landing_rate = [];
                track_start_times = cell(num_mts,1);
                track_dist_to_plus_end = {}; 
                all_landing_dist = []; %not grouped by um, but within 3um
                landing_distance = {};
                all_norm_landing_dist = []; %all_landing_dist but normalized by distance from track of interest to MT end that motor lands towards
                all_mt_landing_dist = []; %all_landing_dist but normalized by total MT length
                all_time_diff_bw_land = []; %stores time between landing events of motors that are within 3um of each other when the later motor lands
                time_diff_bw_land = cell(1,num_mts);
                all_dist_to_minus = [];
                all_dist_to_plus = [];
                pause_durations = {};
                run_durations = {};
                mean_run_vel = {};
                cum_run_durations = [];
                cum_pause_durations = [];
                cum_mean_run_vel = [];
                if zcap == 1
                    landing_dist_by_seg = {};
                    segment_lengths = {};
                    cum_plus_cap_vel = [];
                    cum_plus_gdp_vel = [];
                    cum_seed_vel = [];
                    cum_minus_cap_vel = [];
                    cum_minus_gdp_vel = [];
                    cum_track_start_segment = [];
                    cum_land_rate_by_seg = cell(5,1); %double.empty(5,0); %[];
                    cum_proc_cap_vel = [];
                    cum_proc_gdp_vel = [];
                    cum_proc_seed_vel = [];
                    cum_pause_cap_vel = [];
                    cum_pause_gdp_vel = [];
                    cum_pause_seed_vel = [];
                    cum_mean_proc_cap_vel = [];
                    cum_mean_proc_gdp_vel = [];
                    cum_mean_proc_seed_vel = [];
                end

                %% Analyze
                % profile on
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

                    % FILTER: do not analyze track if it is static (very short run length)
                    if filt_rl == 1 && sqrt((x_tk(end)-x_tk(1))^2+(y_tk(end)-y_tk(1))^2) < min_rl
                        skip_tk = 1;
                        continue
                    end

                    % FILTER: do not analyze track if not near a MT - choose MT along which the given track is, skip track if not near a MT
                    mx_tk = mean(x_tk);
                    my_tk = mean(y_tk);
                    pot_near_mts = zeros(num_mts,1);   
                    for j = 1:num_mts
                        for jj = 1:size(mts{j},1)
                            if ((abs(mts{j}(jj,1)-mx_tk) < 5000) && (abs(mts{j}(jj,2)-my_tk) < 5000)) %initial screen to find close by MTs
                                pot_near_mts(j) = 1; %1 if MT is nearby
                            end
                        end
                    end
                    pot_mt_id = find(pot_near_mts); %indices of potential nearby MTs
                    near_mts = zeros(length(pot_mt_id),1); 
                    for j = 1:length(pot_mt_id)
                        for jj = 1:size(interp_mts{pot_mt_id(j)},1)
                            if ((abs(interp_mts{pot_mt_id(j)}(jj,1)-mx_tk) < mt_dist) && (abs(interp_mts{pot_mt_id(j)}(jj,2)-my_tk) < mt_dist)) %centroid of motor trajectory within mt_dist nm of a position along MT
                                near_mts(j) = 1; %1 if MT is nearby
                            end
                        end
                    end
                    if nnz(near_mts)>0
                        mt_id = pot_mt_id(logical(near_mts)); %indices of nearby MTs
                %         mt_tk = [];
                        if length(mt_id) == 1
                            mt_tk = mt_id;
                            mt_coords = mts{mt_tk};
                        elseif length(mt_id) > 1 %if more than one potential MT nearby
                            mt_error = zeros(nnz(near_mts),1); %vector to store errors
                            for j = 1:nnz(near_mts)
                                [~,uni_ind] = unique(interp_mts{mt_id(j)}(:,1),'stable'); %find repeated x-values in MT
                                %uni_ind = sort(uni_ind);
                                fit_mt = polyfit(interp_mts{mt_id(j)}(uni_ind,1),interp_mts{mt_id(j)}(uni_ind,2),1);%pchip(MTs{chosen_MT}(uni_ind,1),MTs{chosen_MT}(uni_ind,2)); %%%this still "overfits" the MT
                                interp_mt_for_motor = polyval(fit_mt,x_tk);
                                mt_error(j,1) = sum(abs(y_tk - interp_mt_for_motor));
                            end
                            mt_tk = mt_id(mt_error==min(mt_error)); %choose MT that motor trajectory is closest to
                            mt_tk = mt_tk(1); %if more than one MT equally close, choose first MT
                            mt_coords = mts{mt_tk};
                        else
                            skip_tk = 1; %track not near a MT
                            continue
                        end
                    else
                        skip_tk = 1; %track not near a MT
                        continue
                    end

                    %FILTER: do not analyze track if it is along a MT that is to be ignored
                    if any(filt_mt_ind(:) == mt_tk) == 1
                        skip_tk = 1;
                        continue
                    end

                %     % for testing, fit "MT" to trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     pseudo_mt = polyfit(x_tk,y_tk,1);
                %     mt_ycoords = polyval(pseudo_mt,x_tk);
                %     mt_coords = [x_tk, mt_ycoords];
                %     [~,uni_ind] = unique(mt_coords(:,1),'stable'); %find repeated x-values in MT
                %     uni_ind = sort(uni_ind);
                %     mt_coords = mt_coords(uni_ind,:);

                    % If MT is oriented the "wrong" way, this track suggests we reverse MT (actually decided based on tally from all tracks along a given MT)
                    dist_to_start = sqrt((mts{mt_tk}(1,1) - x_tk(1)).^2+(mts{mt_tk}(1,2) - y_tk(1)).^2); %distance from start of MT to start of track
                    dist_to_end = sqrt((mts{mt_tk}(1,1) - x_tk(end)).^2+(mts{mt_tk}(1,2) - y_tk(end)).^2); %distance from start of MT to end of track (if this is smaller, MT is the "wrong" way around
                    if dist_to_start > dist_to_end %motor moves towards start of MT, MT is the "wrong" way
                        flip_mt(mt_tk) = flip_mt(mt_tk)+1; 
                    elseif dist_to_start < dist_to_end 
                        no_flip(mt_tk) = no_flip(mt_tk)+1;
                    end

                    % Find points of track considered to be at end of MT (no distinction made between plus and minus end)
                    motorsq = [x_tk,y_tk];
                    mtendsq = [mt_coords(1,:);mt_coords(end,:)];
                    [points_mtend,dist_mtend] = rangesearch(motorsq,mtendsq,end_dist,'Distance','euclidean'); %points within specified nm of either end of MT (cell with entry for each end) %distance to MT end (cell with entry for each end)
                    if size(points_mtend{1,1},2) ~= 0 || size(points_mtend{2,1},2) ~= 0 %there are points near at least one of the ends
                        ind_mtend = [sort(points_mtend{1,1}), sort(points_mtend{2,1})]; %index of motorsq near MT end
                        loc_mtend = [motorsq(ind_mtend,1), motorsq(ind_mtend,2)]; %x,y values of spots in motorsq near MT end
                        if filt_tk_mt_end == 1
                            cens_tk = 1;
                        else
                            cens_tk = 0;
                        end
                    else
                        ind_mtend = [];
                        loc_mtend = [];
                        cens_tk = 0;
                    end
                    all_ind = 1:1:duration_tk; %all indices (of tk variables)
                    ind_alongmt = setdiff(all_ind,ind_mtend,'stable'); %indices of track (within tk variables) along length of MT

                    %%%%%%%%%%%%%% All filtering is done before this point %%%%%%%%%%%%%%%%
                    ftk = ftk+1; %only advances if track is not filtered
                    censored = [censored; cens_tk]; %cens_tk is 1 if track reaches MT end, 0 if not

                    traj(ftk).mt = mt_tk;
                    traj(ftk).x = x_tk;
                    traj(ftk).y = y_tk;
                    traj(ftk).frames = frame_tk;
                    traj(ftk).duration = duration_tk;
                    traj(ftk).ind_mtend = ind_mtend;
                    traj(ftk).ind_alongmt = ind_alongmt;
                end
                 nfilttracks = ftk;

                % identify plus ends of MTs
                for mttk = 1:num_mts
                    if flip_mt(mttk) > no_flip(mttk)
                        mts{mttk} = flipud(mts{mttk});
                        interp_mts{mttk} = flipud(interp_mts{mttk});
                        %disp(['flipped MT number ',num2str(mttk)])
                    end
                    tot_mt_length = arclength(mts{mttk}(:,1),mts{mttk}(:,2));    
                    mt_lengths{mttk} = tot_mt_length;
                end

                % delete MTs that are to be ignored
                for i = 1:size(filt_mt_ind,2)
                    mts{filt_mt_ind(size(filt_mt_ind,2)-i+1)} = [];
                    interp_mts{filt_mt_ind(size(filt_mt_ind,2)-i+1)} = [];
                end

                motor_on_mt = cell(nfilttracks,1);
                dist_to_bound = {};

                for ftk = 1:nfilttracks
                    mt_tk = traj(ftk).mt;
                    mt_coords = mts{mt_tk};
                    interp_mt_coords = interp_mts{mt_tk};

                    x_tk = traj(ftk).x;
                    y_tk = traj(ftk).y;
                    frame_tk = traj(ftk).frames;
                    duration_tk = traj(ftk).duration;

                    ind_mtend = traj(ftk).ind_mtend;
                    ind_alongmt = traj(ftk).ind_alongmt;

                    %% Find where along MT motor is
                    %find closest points on MT to motor coordinates
                    for j = 1:size(x_tk,1)
                        [nearpoints_motmt,distmotmt_nearpoints] = rangesearch(interp_mts{mt_tk},[x_tk,y_tk],5000,'Distance','euclidean');%,'SortIndices',1); %find nearby points on MT
                        closest_motmtpoint_ind = nearpoints_motmt{j,1}(1);
                        closest_motmtpoint = interp_mts{mt_tk}(closest_motmtpoint_ind,:);
                        motor_on_mt{ftk} = [motor_on_mt{ftk}; closest_motmtpoint];
                    end

                    % start position on MT
                    [~,interp_mt_start_ind] = ismember(motor_on_mt{ftk}(1,:), interp_mts{mt_tk},'rows');
                    interp_mt_start_ind = nonzeros(interp_mt_start_ind);
                    if interp_mt_start_ind == 1
                        dist_to_start = 0;
                    else
                        mt_inds_to_start = 1:1:interp_mt_start_ind;
                        dist_to_start = arclength(interp_mts{mt_tk}(mt_inds_to_start,1),interp_mts{mt_tk}(mt_inds_to_start,2)); %from MT minus-end
                    end

                    tot_mt_length = arclength(mt_coords(:,1),mt_coords(:,2));

                    %% Project track along MT
                    %finding MT/norm(MT)
                    fit_vector = [];
                    motor_vector = [0,0;diff(x_tk),diff(y_tk)];
                    mt_vector = [0,0;diff(mt_coords)];

                    % polynomial fit to vector (gives MT coordinates for each spot localized in track)
                    [~,uni_ind] = unique(interp_mt_coords(:,1),'stable'); %find repeated x-values in MT
                    %uni_ind = sort(uni_ind);
                    fit_mt = polyfit(interp_mt_coords(uni_ind,1),interp_mt_coords(uni_ind,2),1);%pchip(MTs{chosen_MT}(uni_ind,1),MTs{chosen_MT}(uni_ind,2)); %%%this still "overfits" the MT
                    interp_mt_for_motor = polyval(fit_mt,x_tk);
                %     fit_vector(:,1) = x_tk;
                %     fit_vector(:,2) = interp_mt_for_motor;
                    [~,uni_ind2] = unique(mt_coords(:,1),'stable'); %find repeated x-values in MT
                    fit_vector(:,1)=interp1(mt_coords(uni_ind2,1),mt_vector(uni_ind2,1),x_tk,'pchip');
                    fit_vector(:,2)=interp1(mt_coords(uni_ind2,1),mt_vector(uni_ind2,2),x_tk,'pchip');
                    fit_norm = repmat(sqrt(sum(fit_vector.^2,2)),1,2);
                    fit_vector_norm = fit_vector./fit_norm; %normalize

                    % taking dot(track,MT)*MT/norm(MT)
                    delp = zeros(1,numel(x_tk));
                    delp_off = zeros(1,numel(x_tk));
                    for kv = 1:numel(x_tk)
                        delp(kv) = motor_vector(kv,:)*fit_vector_norm(kv,:)'; %displacement projected onto MT axis
                        delp_off(kv) = motor_vector(kv,:)*(fit_vector_norm(kv,:)*[0,1;-1,0])'; %displacement orthogonal to MT axis
                    end
                    position = cumsum(delp); %position of track along MT vector
                    position_off = cumsum(delp_off); %position of track perpendicular to MT vector

                    %% Sliding MSD analysis
                    if length (frame_tk) > l_window + 4
                        tmsd_res = tMSD_2D(x_tk,y_tk,frame_tk,l_window,exp_time,msd_thresh, msd_step,l_min);
                        tmsd_frames = tmsd_res(:,1); %frame of particle (in movie)
                        loc_alpha = tmsd_res(:,2); %local alpha-value
                        proc_frames = tmsd_res(:,3); %was frame marked as processive (1) or paused/diffusive (0)
                        
                        proc_alpha = loc_alpha(proc_frames==1);
                        pause_alpha = loc_alpha(proc_frames==0);
                    end

                    %% Calculate parameters
                    inst_vel = diff(position)./diff(frame_tk.*exp_time);
                    mean_vel = (position(end))/((frame_tk(end)-frame_tk(1)+1)*exp_time);
                    if length (frame_tk) > l_window + 4
                        motorq_mtind = [];
                        run_durations{ftk} = [];
                        pause_durations{ftk} = [];
                        mean_run_vel{ftk} = [];                       
                        proc_vel = []; %inst_vel(proc_frames==1); %diff(position(proc_frames==1))./diff(frame_tk(proc_frames==1).*exp_time);
                        pause_vel = []; %inst_vel(proc_frames==0); %diff(position(proc_frames==0))./diff(frame_tk(proc_frames==0).*exp_time);

                        changes = diff(tmsd_res(:,3)); %-1 means changes processive to paused; 1 means changes paused to processive
                        changeframes =find(changes); %gives the last frame of each segment - these are not absolute frames, but the first frame of a track is always 1
                        segment_ends = [changeframes;size(tmsd_res,1)]; %the last frame of each segment including the last frame - these are not absolute frames, but the first frame of a track is always 1
                        segment_starts = [1; changeframes+1];

                        ind_to_check = sort([segment_ends;segment_starts]);
                        motorqpos = motor_on_mt{ftk};
                        for ik = 1:size(ind_to_check,1)
                            [~,motor_mtind] = ismember(interp_mts{mt_tk}, motorqpos(ind_to_check(ik),:), 'rows');
                            motor_mtind = find(motor_mtind);
                            motorq_mtind = [motorq_mtind;motor_mtind];
                        end
                        
                        pn = 0;
                        rn = 0;
                        if (isempty(changeframes) && ~isempty(tmsd_res)) %motor does not switch between paused and processive
                            if tmsd_res(1,3) == 0
                                pause_durations{ftk} = [pause_durations{ftk}; size(tmsd_res,1)]; %in number of frames
                                pause_vel = diff(position)./diff(frame_tk.*exp_time);
                                pn = pn+1;
                            elseif tmsd_res(1,3) == 1 
                                proc_mt_ind = motorq_mtind(1):1:motorq_mtind(end);
                                if length(proc_mt_ind) > 1
                                    run_durations{ftk} = [run_durations{ftk}; size(tmsd_res,1)]; %in number of frames
                                    proc_disp = arclength(interp_mts{mt_tk}(proc_mt_ind,1),interp_mts{mt_tk}(proc_mt_ind,2));
                                    uprocvel = proc_disp/(size(tmsd_res,1)*exp_time);
                                    mean_run_vel{ftk} = [mean_run_vel{ftk}; uprocvel];
                                    proc_vel = diff(position)./diff(frame_tk.*exp_time);
                                    rn = rn+1;
                                end
                            end
                        elseif(~isempty(changeframes)) %motor does switch
                            spans = [changeframes(1); diff(changeframes); length(changes)+1-changeframes(end)];
                            for sn = 1:(numel(spans))
                                if sn == numel(spans)
                                    if tmsd_res(end,3) == 0
                                        segment_ind = segment_starts(sn):1:size(tmsd_res,1);
                                        pause_vel = [pause_vel, diff(position(segment_ind))./diff(frame_tk(segment_ind)'.*exp_time)];
                                        pause_durations{ftk} = [pause_durations{ftk}; spans(sn)]; %in number of frames
                                        pn = pn  + 1;
                                    elseif tmsd_res(end,3) == 1 
                                        proc_mt_ind = motorq_mtind(2*sn-1):1:motorq_mtind(2*sn);
                                        if length(proc_mt_ind) > 1
                                            segment_ind = segment_starts(sn):1:size(tmsd_res,1);
                                            proc_vel = [proc_vel, diff(position(segment_ind))./diff(frame_tk(segment_ind)'.*exp_time)]; %changed from run_vel
                                            run_durations{ftk} = [run_durations{ftk}; spans(sn)]; %in number of frames
                                            proc_disp = arclength(interp_mts{mt_tk}(proc_mt_ind,1),interp_mts{mt_tk}(proc_mt_ind,2));
                                            uprocvel = proc_disp/(spans(sn)*exp_time);
                                            mean_run_vel{ftk} = [mean_run_vel{ftk}; uprocvel];
                                            rn = rn + 1;
                                        end
                                    end
                                else
                                    if (changes(changeframes(sn)) == 1) %was paused
                                        segment_ind = segment_starts(sn):1:segment_starts(sn+1)-1;
                                        pause_vel = [pause_vel, diff(position(segment_ind))./diff(frame_tk(segment_ind)'.*exp_time)];
                                        %store pause time lengths
                                        pause_durations{ftk} = [pause_durations{ftk}; spans(sn)]; %in number of frames
                                        pn = pn + 1;
                                    elseif (changes(changeframes(sn)) == -1) %was processive
                                        %store run time lengths
                                        proc_mt_ind = motorq_mtind(2*sn-1):1:motorq_mtind(2*sn);
                                        if length(proc_mt_ind) > 1
                                            segment_ind = segment_starts(sn):1:segment_starts(sn+1)-1;
                                            proc_vel = [proc_vel, diff(position(segment_ind))./diff(frame_tk(segment_ind)'.*exp_time)];  %changed from run_vel                                        
                                            run_durations{ftk} = [run_durations{ftk}; spans(sn)]; %in number of frames
                                            proc_disp = arclength(interp_mts{mt_tk}(proc_mt_ind,1),interp_mts{mt_tk}(proc_mt_ind,2));
                                            uprocvel = proc_disp/(spans(sn)*exp_time);
                                            mean_run_vel{ftk} = [mean_run_vel{ftk}; uprocvel];
                                            rn = rn + 1;
                                        end
                                    end
                                end
                            end
                        end
                        num_pauses = pn;
                        num_runs = rn;
                    end 

                    %% Store data
                    traj(ftk).mt = mt_tk;
                    traj(ftk).flip = flip_mt(mt_tk);
                    traj(ftk).x = x_tk;
                    traj(ftk).y = y_tk;
                    traj(ftk).frames = frame_tk;
                    traj(ftk).duration = duration_tk;

                    traj(ftk).ind_alongmt = ind_alongmt;
                    traj(ftk).ind_mtend = ind_mtend;

                    traj(ftk).position = position;
                    traj(ftk).pos_on_interpmt = motor_on_mt{ftk};
                    traj(ftk).start_pos_on_mt = dist_to_start;
                    traj(ftk).mt_coords = mt_coords; %will be flipped if needed for track --> does not necessarily match mts{mt_tk}
                    traj(ftk).interp_mt_coords = interp_mt_coords; %will be flipped if needed for track --> does not necessarily match interp_mts{mt_tk}
                    traj(ftk).offaxis_position = position_off;
                    traj(ftk).mt_length = tot_mt_length;
                    traj(ftk).run_length = position(end);
                    traj(ftk).inst_vel = inst_vel(1,:);

                    if length (frame_tk) > l_window + 4
                        traj(ftk).loc_alpha = loc_alpha;
                        cum_run_durations = [cum_run_durations; run_durations{ftk}];
                        cum_pause_durations = [cum_pause_durations; pause_durations{ftk}];
                        cum_mean_run_vel = [cum_mean_run_vel; mean_run_vel{ftk}];
                        traj(ftk).num_pauses = num_pauses;
                        traj(ftk).num_runs = num_runs;
                    end
                    if length (frame_tk) > l_window + 4
                        traj(ftk).proc_frames = proc_frames;
                    end
                    if length (frame_tk) > l_window + 4 && ~isempty(proc_vel) == 1
                        traj(ftk).proc_vel = proc_vel(1,:);
                    end
                    if length (frame_tk) > l_window + 4 && ~isempty(pause_vel) == 1
                        traj(ftk).pause_vel = pause_vel(1,:);
                    end

                    cum_run_length = [cum_run_length; abs(position(end))]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    cum_censored = censored;
                    cum_mean_vel = [cum_mean_vel; mean_vel];%abs(position(end))/(duration_tk*exp_time)];
                    cum_inst_vel = [cum_inst_vel, inst_vel(1,:)];
                    if length (frame_tk) > l_window + 4
                        cum_loc_alpha = [cum_loc_alpha; loc_alpha];
                    end
                    if length (frame_tk) > l_window + 4 && ~isempty(proc_vel) == 1
                        cum_proc_vel = [cum_proc_vel, proc_vel(1,:)];
                        cum_proc_alpha = [cum_proc_alpha; proc_alpha];
                    end
                    if length (frame_tk) > l_window + 4 && ~isempty(pause_vel) == 1
                        cum_pause_vel = [cum_pause_vel, pause_vel(1,:)];
                        cum_pause_alpha = [cum_pause_alpha; pause_alpha];
                    end
                    cum_mts = [cum_mts; mt_tk];
                    cum_association_time = [cum_association_time; (frame_tk(end)-frame_tk(1)+1)*exp_time];

                end
                % profile viewer

                %% Plot
                if zplot ~= 0
                    for ftk = 1:size(cum_run_length,1)
                        figure(traj_plot), hold on, plot(traj(ftk).x,traj(ftk).y,'.-','Color',cmap(traj(ftk).mt,:)) %plots trajectories
                        figure(kymo_plot), hold on, plot((traj(ftk).frames).*exp_time,traj(ftk).position,'.-','Color',cmap(traj(ftk).mt,:)) %plots "kymographs"
                        %figure(offaxis_plot), hold on, plot((traj(ftk).frames).*exp_time,traj(ftk).offaxis_position,'.-','Color',cmap(ftk,:)) %plots off-axis "kymographs"
                    end
                end    

                %% analyze track start times and position  
                for mttk = 1:num_mts
                    segment_indices = {};
                    ftk_on_mt = find(cum_mts == mttk); %gives indices of cum_mts, which should match that of ftk
                    if ~isempty(mts{mttk}) == 1
                        tot_mt_length = arclength(mts{mttk}(:,1),mts{mttk}(:,2));
                    end
                    
                    landing_rate = [landing_rate; numel(ftk_on_mt)/((tot_mt_length/1000)*(num_frames*exp_time)*motor_concentration)]; %[um-1s-1nM-1]
                        
                    landing_distance{mttk} = cell(max(ftk_on_mt),max(ftk_on_mt)); %cell(size(ftk_on_mt,1),size(ftk_on_mt,1));
                    landing_dist_by_seg{mttk} = cell(5,1);
                    if ~isempty(ftk_on_mt)
                        if zplot ~= 0
                            figure, hold on
                        end
                        
                        for j=1:length(ftk_on_mt)      
                            pos_on_mt = traj(ftk_on_mt(j)).position+traj(ftk_on_mt(j)).start_pos_on_mt;

                            if zplot ~= 0
                                %plot((traj(ftk_on_mt(j)).frames).*exp_time,pos_on_mt,'.-','Color',cmap(mttk,:))%ftk_on_mt(j),:)) %plots "kymograph" for each MT
                                plot((traj(ftk_on_mt(j)).frames).*exp_time,pos_on_mt,'.-','Color',cmap(mttk,:))%ftk_on_mt(j),:)) %plots "kymograph" for each MT
                                %plot((traj(ftk_on_mt(j)).frames).*exp_time,traj(ftk_on_mt(j)).position,'.-','Color',cmap(ftk_on_mt(j),:)) %plots "kymographs"
                            end

                           % start position on MT
                            [~,interp_mt_end_ind] = ismember(traj(ftk_on_mt(j)).pos_on_interpmt(1,:),interp_mts{mttk}, 'rows');
                            interp_mt_end_ind = nonzeros(interp_mt_end_ind);
                            if interp_mt_end_ind == size(interp_mts{mttk},1)
                                dist_to_end = 0;
                            else
                                mt_inds_to_end = interp_mt_end_ind:1:size(interp_mts{mttk},1); % indices to plus-end
                                dist_to_end = arclength(interp_mts{mttk}(mt_inds_to_end,1),interp_mts{mttk}(mt_inds_to_end,2));
                            end
                            landing_dist_to_mt_end = dist_to_end; %distance from start of track to plus-end of MT
                            normalized_landing_pos = traj(ftk_on_mt(j)).start_pos_on_mt/tot_mt_length;

                            %save position on MT and landing info
                            traj(ftk_on_mt(j)).position_on_mt = pos_on_mt;
                            traj(ftk_on_mt(j)).normalized_pos_on_mt = pos_on_mt./tot_mt_length;
                            cum_norm_landing_pos = [cum_norm_landing_pos; normalized_landing_pos];
                            cum_landing_dist_to_mt_end = [cum_landing_dist_to_mt_end; landing_dist_to_mt_end];

                            %landing rate with time
                            track_start_times{mttk} = [track_start_times{mttk}; traj(ftk_on_mt(j)).frames(1)*exp_time]; %[s]
                            
                            %% analyzing proximity of landing of other motors to the track in SPACE
                            frames_to_check = traj(ftk_on_mt(j)).frames;
                            motorqpos = traj(ftk_on_mt(j)).pos_on_interpmt;
                            other_traj_ind = setdiff(ftk_on_mt,ftk_on_mt(j));

%                             if tot_mt_length >= 6000 && traj(ftk_on_mt(j)).position_on_mt(1) <= 6000 %lands within first 6um of a MT at least 6um long
%                                 for jk = 1:size(other_traj_ind,1)
%                                     if traj(other_traj_ind(jk)).position_on_mt <= 6000 %also lands within first 6um of MT
%                                         all_time_diff_bw_land = [all_time_diff_bw_land; (traj(other_traj_ind(jk)).frames(1)-traj(ftk_on_mt(j)).frames(1))*exp_time]; %[s]
%                                     end
%                                 end
%                             end
                            
                            for framek = 1:size(frames_to_check,1)
                                [~,motorq_mtind] = ismember(interp_mts{mttk}, motorqpos(framek,:), 'rows');
                                motorq_mtind = find(motorq_mtind);
                                if motorq_mtind == 1
                                    motor_to_minus = 10;
                                else
                                    motor_to_minus_ind = 1:1:motorq_mtind;
                                    motor_to_minus = arclength(interp_mts{mttk}(motor_to_minus_ind,1),interp_mts{mttk}(motor_to_minus_ind,2)); %distance to minus end
                                end
                                if motorq_mtind == size(interp_mts{mttk},1)
                                    motor_to_plus = 10;
                                else
                                    motor_to_plus_ind = motorq_mtind:1:size(interp_mts{mttk},1);
                                    motor_to_plus = arclength(interp_mts{mttk}(motor_to_plus_ind,1),interp_mts{mttk}(motor_to_plus_ind,2)); %distance to plus end
                                end
                                if motor_to_minus >= 3000 && motor_to_plus >= 3000
                                    for jk = 1:size(other_traj_ind,1)
                                        if ~isempty(other_traj_ind) && traj(other_traj_ind(jk)).frames(1) == frames_to_check(framek)
                                            [~,land_mtind] = ismember(interp_mts{mttk}, traj(other_traj_ind(jk)).pos_on_interpmt(1,:), 'rows');
                                            land_mtind = find(land_mtind);
                                            if land_mtind == motorq_mtind
                                                dist_bw = 10; %if they are assigned to the same interp_mt index, they are at most 10nm apart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRY THIS WITH 0 
                                                neg_dist = 0; %assume towards plus end
                                                norm_dist = motor_to_plus;
                                            else
                                                ind_temp = [land_mtind,motorq_mtind];
                                                ind_temp = sort(ind_temp);
                                                if ind_temp(1) == land_mtind %motor landed behind motorq (i.e. between motor and minus-end)
                                                    neg_dist = 1;
                                                    norm_dist = motor_to_minus;
                                                else %motor landed ahead of motorq (i.e. between motor and plus-end)
                                                    neg_dist = 0;
                                                    norm_dist = motor_to_plus;
                                                end
                                                dist_ind = ind_temp(1):1:ind_temp(end);
                                                dist_bw = arclength(interp_mts{mttk}(dist_ind,1),interp_mts{mttk}(dist_ind,2));
                                                if neg_dist == 1
                                                    dist_bw = -1*dist_bw;
                                                end
                                            end
                                            norm_dist_bw = dist_bw/norm_dist;
                                            mt_dist_bw = dist_bw/tot_mt_length;
                                            all_norm_landing_dist = [all_norm_landing_dist;norm_dist_bw];
                                            all_mt_landing_dist = [all_mt_landing_dist;mt_dist_bw];
%                                            all_time_diff_bw_land = [all_time_diff_bw_land; (frames_to_check(framek)-frames_to_check(1))*exp_time]; %[s]
                                            if abs(dist_bw) <= 3000 %lands within 3um of existing track
                                                all_landing_dist = [all_landing_dist; dist_bw]; %within 3um
                                                landing_distance{mttk}{ftk_on_mt(j),other_traj_ind(jk)} = [frames_to_check(framek), dist_bw]; %{mttk}{ftk_on_mt(j),other_traj_ind(jk)}
                                          end
                                        end
                                    end
                                end
                                all_dist_to_plus = [all_dist_to_plus; motor_to_plus];
                                all_dist_to_minus = [all_dist_to_minus; motor_to_minus];
                            end
                           
                %                 figure,trackstarts=gcf; %initialize figure
                %                 [tkstart_n, tkstart_edges]=histcounts(cum_run_length, 'BinWidth', tkstart_binwidth, 'Normalization', 'pdf');
                %                 nhist_tkstart=tkstart_n;
                %                 xhist_tkstart=tkstart_edges+(tkstart_binwidth/2);
                %                 xhist_tkstart(end)=[];
                %                 figure(trackstarts), hold on
                %                 bar(xhist_tkstart,nhist_tkstart)
                %                 xlabel('Time of landing (s)'), ylabel('Probability density'), title([motor,' ', mt_type,' Track start times for MT number ', num2str(mttk)])
                        end
                        
                        %% analyzing proximity of landing of other motors to the track in TIME
                        if tot_mt_length >= 6000 %MT at least 6um long
                            traj_in_mt_sec = [];
                            for i=1:length(ftk_on_mt)
                                if traj(ftk_on_mt(i)).position_on_mt(1) <= 6000 %find indices of all tracks that start within first 6um of MT
                                    traj_in_mt_sec = [traj_in_mt_sec; ftk_on_mt(i)];
                                end
                            end
                            track_start_t_on_mt = [];
                            for i=1:length(traj_in_mt_sec)
                                track_start_t_on_mt(i) = traj(traj_in_mt_sec(i)).frames(1); %corresponding start times of tracks
                            end
                            time_diff_bw_land{mttk} = diff(track_start_t_on_mt)'.*exp_time; %[s]
                            all_time_diff_bw_land = [all_time_diff_bw_land; diff(track_start_t_on_mt)'.*exp_time]; %[s]
                            %PUT MATCHING DISTANCE B/W LANDING
                            %HERE!!! %%%%%%%%%
                        end
        
                        if zplot ~= 0
                            xlabel('time (s)'), ylabel('position along MT (nm)'), title(['Kymographs for MT number ', num2str(mttk), ' (MT length ', num2str(tot_mt_length), 'nm)'])
                            xlim([0 num_frames*exp_time])%, ylim([0 mt_length])
                            hold off
                        end

                        %% Identify which parts of track are on GMP-CPP parts of MT and which are on GDP parts
                        if zcap == 1
                            %order segment boundaries based on distance from plus-end of MT
                            if ~isempty(boundaries_on_mt{mttk})
                                boundaries_on_mt_old = boundaries_on_mt;
                                for j=1:size(boundaries_on_mt{mttk},1)
                                    [~,interp_mt_bound_ind] = ismember(boundaries_on_mt{mttk}(j,:),interp_mts{mttk}, 'rows');
                                    interp_mt_bound_ind = nonzeros(interp_mt_bound_ind);
                                    mt_inds_to_bound = interp_mt_bound_ind:1:size(interp_mts{mttk},1);
                                    dist_to_bound{mttk}(j,1) = arclength(interp_mts{mttk}(mt_inds_to_bound,1),interp_mts{mttk}(mt_inds_to_bound,2)); %distance to plus-end
                                    dist_to_bound{mttk}(j,2) = j;
                                end
                                dist_to_bound{mttk} = sortrows(dist_to_bound{mttk},1);
                                for j=1:size(boundaries_on_mt{mttk},1)
                                    boundaries_on_mt{mttk}(j,:) = boundaries_on_mt_old{mttk}(dist_to_bound{mttk}(j,2),:); %boundaries ordered from closest to furthest from plus-end of MT 
                                end
                            end
                            dist_along_mt = [0;dist_to_bound{mttk}(:,1);tot_mt_length];

                            % find length of each segment
                            segment_lengths{mttk} = abs(diff(dist_along_mt)); %either 3 or 5 entries in order +cap, +gdp, seed(, -gdp, -cap)

                            if size(dist_along_mt,1) == 4 
                                segment_annotate = 3;
                            elseif size(dist_along_mt,1) == 6
                                segment_annotate = 5;
                            else
                                segment_annotate = 0;
                                disp(['MT number ', num2str(mttk), ' does not have 3 or 5 segments. Retrace MT and/or segment boundaries']);
                            end

                            %2 cases, either we have 3 (seed, GDP, cap) or 5 (-cap, GDP, seed, GDP, +cap) segments
                            if segment_annotate ~= 0
                                for i=1:length(ftk_on_mt) %each track on this MT
                                    %find distance of each point in motor trajectory to MT minus-end
                                    track_on_interp_mt = traj(ftk_on_mt(i)).pos_on_interpmt;
                                    for k = 1:size(track_on_interp_mt,1)
                                        [~,interp_mt_point_ind] = ismember(track_on_interp_mt(k,:),interp_mts{mttk}, 'rows');
                                        interp_mt_point_ind = nonzeros(interp_mt_point_ind);
                                        if interp_mt_point_ind == size(interp_mts{mttk},1)
                                            track_dist_to_plus_end{mttk,ftk_on_mt(i)}(k,1) = 0;
                                        else
                                            mt_inds_to_point = interp_mt_point_ind:1:size(interp_mts{mttk},1);
                                            track_dist_to_plus_end{mttk,ftk_on_mt(i)}(k,1) = arclength(interp_mts{mttk}(mt_inds_to_point,1),interp_mts{mttk}(mt_inds_to_point,2));
                                        end
                                    end

                                    %see which segment of MT each point in track is on
                                    for k = 1:segment_annotate
                                        [ind,~] = find(track_dist_to_plus_end{mttk,ftk_on_mt(i)}(:,1)>dist_along_mt(k) & track_dist_to_plus_end{mttk,ftk_on_mt(i)}(:,1)<dist_along_mt(k+1));
                                        if k == 1 %is on +cap
                                            plus_cap_ind = ind;
                                        elseif k == 2 %is on gdp
                                            plus_gdp_ind = ind;
                                        elseif k == 3 %is on seed
                                            seed_ind = ind;
                                        elseif k == 4 %is on gdp
                                            minus_gdp_ind = ind;
                                        elseif k == 5 %is on -cap
                                            minus_cap_ind = ind;
                                        end
                                        segment_indices{mttk}{ftk_on_mt(i),k} = ind;
                                        if ~isempty(ind == 1)
                                            track_start_segment = k;
                                        end
                                    end
                                                                        
                                    %analyze velocities on different segments
                                    plus_cap_vel = diff(traj(ftk_on_mt(i)).position(segment_indices{mttk}{ftk_on_mt(i),1}))./diff(traj(ftk_on_mt(i)).frames(segment_indices{mttk}{ftk_on_mt(i),1}).*exp_time);
                                    if ~isempty(plus_cap_vel)
                                        plus_cap_vel = plus_cap_vel(1,:);
                                    end
                                    plus_gdp_vel = diff(traj(ftk_on_mt(i)).position(segment_indices{mttk}{ftk_on_mt(i),2}))./diff(traj(ftk_on_mt(i)).frames(segment_indices{mttk}{ftk_on_mt(i),2}).*exp_time);
                                    if ~isempty(plus_gdp_vel)
                                        plus_gdp_vel = plus_gdp_vel(1,:);
                                    end
                                    seed_vel = diff(traj(ftk_on_mt(i)).position(segment_indices{mttk}{ftk_on_mt(i),3}))./diff(traj(ftk_on_mt(i)).frames(segment_indices{mttk}{ftk_on_mt(i),3}).*exp_time);
                                    if ~isempty(seed_vel)
                                        seed_vel = seed_vel(1,:);
                                    end
                                    if segment_annotate == 5
                                        minus_cap_vel = diff(traj(ftk_on_mt(i)).position(segment_indices{mttk}{ftk_on_mt(i),5}))./diff(traj(ftk_on_mt(i)).frames(segment_indices{mttk}{ftk_on_mt(i),5}).*exp_time);
                                        if ~isempty(minus_cap_vel)
                                            minus_cap_vel = minus_cap_vel(1,:);
                                        end
                                        minus_gdp_vel = diff(traj(ftk_on_mt(i)).position(segment_indices{mttk}{ftk_on_mt(i),4}))./diff(traj(ftk_on_mt(i)).frames(segment_indices{mttk}{ftk_on_mt(i),4}).*exp_time);
                                        if ~isempty(minus_gdp_vel)
                                            minus_gdp_vel = minus_gdp_vel(1,:);
                                        end
                                    else
                                        minus_cap_vel = [];
                                        minus_gdp_vel = [];
                                    end
                                    cum_plus_cap_vel = [cum_plus_cap_vel, plus_cap_vel];
                                    cum_plus_gdp_vel = [cum_plus_gdp_vel, plus_gdp_vel];
                                    cum_seed_vel = [cum_seed_vel, seed_vel];
                                    cum_minus_cap_vel = [cum_minus_cap_vel, minus_cap_vel];
                                    cum_minus_gdp_vel = [cum_minus_gdp_vel, minus_gdp_vel];
                                    cum_track_start_segment = [cum_track_start_segment, track_start_segment];

                                    %%%%%%
                                    if zplot ~= 0
                                        figure
                                        plot(mts{mttk}(:,1),mts{mttk}(:,2),'-','LineWidth',3)
                                        hold on
                                        plot(traj(ftk_on_mt(i)).x,traj(ftk_on_mt(i)).y,'.-','LineWidth',3,'MarkerSize',20)
                                        plot(track_on_interp_mt(:,1),track_on_interp_mt(:,2),'r.-','LineWidth',3,'MarkerSize',20)
                                        plot(track_on_interp_mt(plus_cap_ind,1),track_on_interp_mt(plus_cap_ind,2),'m.-','LineWidth',3,'MarkerSize',20)
                                        plot(track_on_interp_mt(plus_gdp_ind,1),track_on_interp_mt(plus_gdp_ind,2),'g.-','LineWidth',3,'MarkerSize',20)
                                        plot(track_on_interp_mt(seed_ind,1),track_on_interp_mt(seed_ind,2),'k.-','LineWidth',3,'MarkerSize',20)
                                        scatter(mts{mttk}(:,1),mts{mttk}(:,2),28,'k', 'filled')
                                        plot(mts{mttk}(1,1),mts{mttk}(1,2),'k*','MarkerSize',20)
                                        plot(track_on_interp_mt(1,1),track_on_interp_mt(1,2),'k*','MarkerSize',20)
                                        if zcap ==1
                                            scatter(boundaries_on_mt{mttk}(:,1),boundaries_on_mt{mttk}(:,2),28,'k','filled')
                                            scatter(boundaries_on_mt{mttk}(1,1),boundaries_on_mt{mttk}(1,2),28,'m','filled')
                                            scatter(boundaries_on_mt{mttk}(2,1),boundaries_on_mt{mttk}(2,2),28,'g','filled')
                                        end
                                    end
                                    %save results
                                    traj(ftk_on_mt(i)).plus_cap_vel = plus_cap_vel;
                                    traj(ftk_on_mt(i)).plus_gdp_vel = plus_gdp_vel;
                                    traj(ftk_on_mt(i)).seed_vel = seed_vel;
                                    traj(ftk_on_mt(i)).minus_gdp_vel = minus_gdp_vel;
                                    traj(ftk_on_mt(i)).minus_cap_cel = minus_cap_vel;
                                    traj(ftk_on_mt(i)).track_start_segment = track_start_segment;  

                                    tot_ind = numel(traj(ftk_on_mt(i)).frames);

                                    traj(ftk_on_mt(i)).plus_cap_ind = zeros(tot_ind,1);
                                    traj(ftk_on_mt(i)).plus_gdp_ind = zeros(tot_ind,1);
                                    traj(ftk_on_mt(i)).seed_ind = zeros(tot_ind,1);
                                    
                                    %if segment_annotate == 5
                                    %if ~isempty(minus_gdp_vel)
                                        traj(ftk_on_mt(i)).minus_gdp_ind = zeros(tot_ind,1);
%                                     else
%                                         traj(ftk_on_mt(i)).minus_gdp_ind = [];
                                    %end
                                    %if ~isempty(minus_cap_vel)
                                        traj(ftk_on_mt(i)).minus_cap_ind = zeros(tot_ind,1);
%                                     else
%                                         traj(ftk_on_mt(i)).minus_cap_ind = [];
                                    %end
                                    %end

                                    traj(ftk_on_mt(i)).plus_cap_ind(plus_cap_ind) = 1;
                                    traj(ftk_on_mt(i)).plus_gdp_ind(plus_gdp_ind) = 1;
                                    traj(ftk_on_mt(i)).seed_ind(seed_ind) = 1;
                                    if ~isempty(minus_gdp_vel)
                                        traj(ftk_on_mt(i)).minus_gdp_ind(minus_gdp_ind) = 1;
                                    end
                                    if ~isempty(minus_cap_vel)
                                        traj(ftk_on_mt(i)).minus_cap_ind(minus_cap_ind) = 1;
                                    end
                                    
                                    
                                    %%% ADD STUFF HERE ABOUT PROC AND GDP
                                    %%% ETC.
                                    if length(traj(ftk_on_mt(i)).frames) > l_window + 4
                                        temp_a = isempty(traj(ftk_on_mt(i)).proc_frames); %exist('traj(ftk_on_mt(i)).proc_frames');
                                        if temp_a ~= 1 %0
                                        %if ~isempty(traj(ftk_on_mt(i)).proc_frames) == 1 %length (traj(ftk_on_mt(i)).frames) > l_window + 4
                                        pause_ind = ones(size(traj(ftk_on_mt(i)).proc_frames,1),1);
                                        pause_ind = pause_ind -traj(ftk_on_mt(i)).proc_frames;
                                            %if ~isempty(minus_cap_vel) == 1
                                            if segment_annotate == 5
                                                proc_cap_ind = [traj(ftk_on_mt(i)).plus_cap_ind+traj(ftk_on_mt(i)).minus_cap_ind].* traj(ftk_on_mt(i)).proc_frames;
                                                proc_gdp_ind = [traj(ftk_on_mt(i)).plus_gdp_ind+traj(ftk_on_mt(i)).minus_gdp_ind].* traj(ftk_on_mt(i)).proc_frames;
                                                proc_seed_ind = traj(ftk_on_mt(i)).seed_ind.* traj(ftk_on_mt(i)).proc_frames;
                                                pause_cap_ind = [traj(ftk_on_mt(i)).plus_cap_ind+traj(ftk_on_mt(i)).minus_cap_ind].* pause_ind;
                                                pause_gdp_ind = [traj(ftk_on_mt(i)).plus_gdp_ind+traj(ftk_on_mt(i)).minus_gdp_ind].* pause_ind;
                                                pause_seed_ind = traj(ftk_on_mt(i)).seed_ind.* pause_ind;
                                            else
                                                proc_cap_ind = [traj(ftk_on_mt(i)).plus_cap_ind].* traj(ftk_on_mt(i)).proc_frames;
                                                proc_gdp_ind = [traj(ftk_on_mt(i)).plus_gdp_ind].* traj(ftk_on_mt(i)).proc_frames;
                                                proc_seed_ind = traj(ftk_on_mt(i)).seed_ind.* traj(ftk_on_mt(i)).proc_frames;
                                                pause_cap_ind = [traj(ftk_on_mt(i)).plus_cap_ind].* pause_ind;
                                                pause_gdp_ind = [traj(ftk_on_mt(i)).plus_gdp_ind].* pause_ind;
                                                pause_seed_ind = traj(ftk_on_mt(i)).seed_ind.* pause_ind;
                                            end
                                            traj(ftk_on_mt(i)).proc_cap_ind = proc_cap_ind;
                                            traj(ftk_on_mt(i)).proc_gdp_ind = proc_gdp_ind;
                                            traj(ftk_on_mt(i)).proc_seed_ind = proc_seed_ind;
                                            traj(ftk_on_mt(i)).pause_cap_ind = pause_cap_ind;
                                            traj(ftk_on_mt(i)).pause_gdp_ind = pause_gdp_ind;
                                            traj(ftk_on_mt(i)).pause_seed_ind = pause_seed_ind;

                                            % look at processive, paused velocities in different segments
                                            % identify first frames of discontinuous segments
                                            % (if motor switches from one condition and back to it there will be >1 changept, with the index of the first point of the new segment identified in ipt

                                            seg_type = zeros(tot_ind,1);
                                            seg_type(proc_cap_ind==1) = 1;
                                            seg_type(proc_gdp_ind==1) = 2;
                                            seg_type(proc_seed_ind==1) = 3;
                                            seg_type(pause_cap_ind==1) = 4;
                                            seg_type(pause_gdp_ind==1) = 5;
                                            seg_type(pause_seed_ind==1) = 6;

                                            seg_type_vel = cell(6,1);
                                            for lmn = 1:6
                                                seg_type_vel{lmn,1} = diff(traj(ftk_on_mt(i)).position(seg_type==lmn))'./(diff(traj(ftk_on_mt(i)).frames(seg_type==lmn)).*exp_time);
                                            end
                                            
                                            %store data
                                            cum_proc_cap_vel = [cum_proc_cap_vel;seg_type_vel{1,1}];
                                            cum_proc_gdp_vel = [cum_proc_gdp_vel;seg_type_vel{2,1}];
                                            cum_proc_seed_vel = [cum_proc_seed_vel;seg_type_vel{3,1}];
                                            cum_pause_cap_vel = [cum_pause_cap_vel;seg_type_vel{4,1}];
                                            cum_pause_gdp_vel = [cum_pause_gdp_vel;seg_type_vel{5,1}];
                                            cum_pause_seed_vel = [cum_pause_seed_vel;seg_type_vel{6,1}];

                                            %find mean velocities of processive segments on different lattices
                                            [mu_proccapvels] = mean_segment_vel(motor_on_mt{ftk_on_mt(i)}, traj(ftk_on_mt(i)).position, traj(ftk_on_mt(i)).frames, exp_time, interp_mts{mttk}, proc_cap_ind);
                                            [mu_procgdpvels] = mean_segment_vel(motor_on_mt{ftk_on_mt(i)}, traj(ftk_on_mt(i)).position, traj(ftk_on_mt(i)).frames, exp_time, interp_mts{mttk}, proc_gdp_ind);
                                            [mu_procseedvels] = mean_segment_vel(motor_on_mt{ftk_on_mt(i)}, traj(ftk_on_mt(i)).position, traj(ftk_on_mt(i)).frames, exp_time, interp_mts{mttk}, proc_seed_ind);
                                            
                                            %store data
                                            cum_mean_proc_cap_vel = [cum_mean_proc_cap_vel;mu_proccapvels];
                                            cum_mean_proc_gdp_vel = [cum_mean_proc_gdp_vel;mu_procgdpvels];
                                            cum_mean_proc_seed_vel = [cum_mean_proc_seed_vel;mu_procseedvels];
                                        end
                                    end
                                end
                                
                                
                                
                                
                                
                                %analyze landing distance to other motors based on segments
                                % landing_distance{mttk}{ftk_on_mt(j),other_traj_ind(jk)} = [framek, dist_bw];
                                %segment_indices{mttk}{ftk_on_mt(i),k} with
                                %k = 1 plus cap, k=2 plus gdp, k=3 seed,
                                %...
                                if ~isempty(landing_distance{mttk}) == 1
                                    for traji = 1:size(ftk_on_mt,1) %each track on MT
                                        other_traj_ind = setdiff(ftk_on_mt,ftk_on_mt(traji));
                                        for trajj = size(other_traj_ind, 1) %each "other" track on MT
                                            if ~isempty(other_traj_ind) == 1 && ~isempty(landing_distance{mttk}{ftk_on_mt(traji),other_traj_ind(trajj)}) == 1
                                                frameq = landing_distance{mttk}{ftk_on_mt(traji),other_traj_ind(trajj)}(1,1);
                                                for segk = 1:size(segment_indices{mttk},2) %each segment (1:3 or 1:5)
                                                    if ismember(frameq,traj(ftk_on_mt(traji)).frames(segment_indices{mttk}{ftk_on_mt(traji),segk})) == 1
                                                        motqseg = segk;
                                                        checkseg = 1;
                                                    else
                                                        checkseg = 0;
                                                        motqseg = -1;
                                                    end
                                                    if ismember(frameq,traj(other_traj_ind(trajj)).frames(segment_indices{mttk}{other_traj_ind(trajj),segk})) == 1
                                                        otherseg = segk;
                                                        checkseg = 1;
                                                    else
                                                        checkseg = 0;
                                                        otherseg = -2;
                                                    end
                                                    if checkseg == 1 && motqseg == otherseg %motors present on same seg %%%%% currently not checking whether segment itself is long enough!!
                                                        landing_dist_by_seg{mttk}{motqseg} = [landing_dist_by_seg{mttk}{motqseg},landing_distance{mttk}{ftk_on_mt(traji),other_traj_ind(trajj)}(1,2)]; %{mttk}{ftk_on_mt(i),other_traj_ind(j)}
                                                    end
                                                end 
                                            end
                                        end
                                    end
                                end
                            end
                            land_rate_by_seg = double.empty(5,0);
                            mtedges = 0.5:1:(0.5+segment_annotate);
                            [land_num_by_seg, ~]=histcounts(cum_track_start_segment, 'BinEdges', mtedges, 'Normalization', 'count');
                            land_rate_by_seg = land_num_by_seg'./((segment_lengths{mttk}./1000).*num_frames.*exp_time.*motor_concentration);
                            for segi = 1:segment_annotate
                                cum_land_rate_by_seg{segi} = [cum_land_rate_by_seg{segi}, land_rate_by_seg(segi)];
                            end
                        end
                    end
                end
                
%                 temp_a = exist('land_rate_by_seg');
%                 if zcap == 1 && temp_a == 1
%                     cum_land_rate_by_seg = [cum_land_rate_by_seg, land_rate_by_seg];
%                 end
                
                %% Plot parameters to check data
                if zplot ~=0
                %     figure,totrl=gcf; %initialize figure
                %     [trl_n, trl_edges]=histcounts(cum_run_length, 'BinWidth', rl_binwidth, 'Normalization', 'pdf');
                %     nhist_trl=trl_n;
                %     xhist_trl=trl_edges+(rl_binwidth/2);
                %     xhist_trl(end)=[]; 
                %     figure(totrl), hold on 
                %     bar(xhist_trl,nhist_trl)
                %     xlabel('Total run length (nm)'), ylabel('Probability density'), title([motor,' ', mt_type,' Total run length histogram'])
                %     opt = statset('mlecustom');
                %     opt = statset(opt,'FunValCheck','off','MaxIter',1e5,'MaxFunEvals',1e5,'Display','iter','TolFun',10e-20);
                %     p0 = [1000];
                %     loL = [300];
                %     upL = [3000];
                %     [estimRL, pciRL] = mle(cum_run_length,'Distribution','exponential','Censoring',censored,'start',p0,'Options',opt,'LowerBound', loL);%'UpperBound', upL) %
                %     %plot fit
                %     yRL = exppdf(xhist_trl, estimRL);
                %     yRLlo = exppdf(xhist_trl, pciRL(1));
                %     yRLhi = exppdf(xhist_trl, pciRL(2));
                %     plot(xhist_trl,0.8.*yRL,'r-');
                %     plot(xhist_trl,0.8.*yRLlo,'r.');
                %     plot(xhist_trl,0.8.*yRLhi,'r.');
                %     hold off

                end

                %% Save data
                if zsave ~= 0
                    %save_dirname =strcat('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data\results'); %windows
                    save_dirname =strcat('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/in_vitro_data/results'); %mac
                    save_filename = ['post_particle_tracking','_',date,'_',motor,'_',mt_type,'_',num2str(filenum),'.mat'];
                    if zcap == 1
                        save(fullfile(save_dirname,save_filename),'mts','interp_mts','traj','track_start_times','cum_run_length','cum_censored', 'cum_mean_vel','cum_inst_vel','cum_proc_vel','cum_pause_vel','cum_loc_alpha','cum_proc_alpha','cum_pause_alpha','cum_association_time', 'cum_norm_landing_pos', 'cum_landing_dist_to_mt_end','boundaries_on_mt','segment_lengths','cum_plus_cap_vel','cum_plus_gdp_vel', 'cum_seed_vel', 'cum_minus_cap_vel','cum_minus_gdp_vel','cum_track_start_segment','mt_lengths','all_landing_dist','all_norm_landing_dist','all_mt_landing_dist','all_dist_to_plus','all_dist_to_minus','all_time_diff_bw_land','time_diff_bw_land','landing_dist_by_seg','segment_indices','landing_rate', 'run_durations', 'pause_durations', 'mean_run_vel', 'cum_run_durations','cum_pause_durations','cum_mean_run_vel','cum_land_rate_by_seg','cum_proc_cap_vel','cum_proc_gdp_vel','cum_proc_seed_vel','cum_pause_cap_vel','cum_pause_gdp_vel','cum_pause_seed_vel','cum_mean_proc_cap_vel','cum_mean_proc_gdp_vel','cum_mean_proc_seed_vel','motor_concentration')%,'all_land_dist','tot_xhist_landdist','tot_nhist_landdist','x_ld_all','n_ld_all','xhist_ld', 'nhist_ld')
                        %writematrix(proc_vel,[save_dirname,save_filename,'proc_vel','.csv'])
                        %save_land_rates = landing_rates(find(~cellfun('isempty', mts)));
                        %writematrix(save_land_rates,[save_dirname,save_filename,'land_rates','.csv'])
                    else
                        save(fullfile(save_dirname,save_filename),'mts','interp_mts','traj','track_start_times','cum_run_length','cum_censored', 'cum_mean_vel','cum_inst_vel','cum_proc_vel','cum_pause_vel','cum_loc_alpha','cum_proc_alpha','cum_pause_alpha','cum_association_time', 'cum_norm_landing_pos', 'cum_landing_dist_to_mt_end','mt_lengths','all_landing_dist','all_norm_landing_dist','all_mt_landing_dist','all_dist_to_plus','all_dist_to_minus','all_time_diff_bw_land','time_diff_bw_land','landing_rate', 'run_durations', 'pause_durations', 'mean_run_vel', 'cum_run_durations','cum_pause_durations','cum_mean_run_vel','motor_concentration')%,'all_land_dist','tot_xhist_landdist','tot_nhist_landdist','x_ld_all','n_ld_all','xhist_ld', 'nhist_ld')
                        %writematrix(proc_vel,[save_dirname,save_filename,'proc_vel','.csv'])
                        %save_land_rates = landing_rates(find(~cellfun('isempty', mts)));
                        %writematrix(save_land_rates,[save_dirname,save_filename,'land_rates','.csv'])
                    end
                end
            end
        end
    end
end