
%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski and C Chen 2019-11-05
% Projection along MT taken from AG Hendricks Lab (McGill)

%%
% TO DO:
% in filter_mts, filter out MTs that run out of edge of FOV
% sliding MSD: window size and 2D
% landing rate per MT with time (also for whole MT for whole movie div by L MT and L mov)

%% This code reads in output from the DoM Utrecht plug-in for detection and tracking of particles for a single movie.
% It reads in all the data and re-organizes it into a format used by cumulative_track_analysis_2.
% It can be used to visualize information about the trajectories and trial parameters.
% Filtering done in this code can be controlled using Options and Parameters

% To use this file you will need the following files:
% Results from DoM Utrecht particle localization and tracking (csv; file name must start with DoM)
% MT x,y positions (csv; file name must start with MT)
%

clear all, close all
%% Options (make 0 to NOT perform related action, 1 to perform)
zplot = 1;
zsave = 0;
% Filtering
filt_cross_mt = 1; %ignore any tracks on MTs that are too close to another MT

set(0,'DefaultFigureWindowStyle','docked')
%addpath(' ')

%% Parameters
% From imaging:
pixel_size = 64.0; %pixel size of camera [nm]
exp_time = 0.1; %exposure time [s]
num_frames = 600; %number of frames in movie [frames]
num_pix_x = 512;
num_pix_y = 512;

% For analysis:
min_duration = 5; %minimum length of track to be analyzed [frames]
end_dist = 200; %maximum distance to first/last point of a MT for a spot localization to be considered to be at MT end [nm]
mt_dist = 400; %maximum distance to some MT for track to considered on a MT and thus analyzed [nm]
analyze_mt_num = -1; %specify which MT to analyze, based on MT id in text file (these start at 0); if this is -1, all MTs will be analyzed
mt_cross_dist = 200; %distance between points on two MTs for them to be considered too close/crossing --> filtered out [nm]
l_window = 7; %number of frames to average for sliding MSD analysis
msd_thresh = 1.1; %alpha-value above which is processive, below which is paused
msd_step = 0.3; %minimum threshold for findchangepts function; minimum improvement in residual error; changepoints currently based on mean values, but can use rms, std, mean and slope
l_min = 3; %minimum distance between two changepoints - smallest duration of pause/run

% For plotting:
rl_binwidth = 100; %bin width for run length histograms

%% Movie to analyze
motor = 'kif1a';
mt_type = '1cycle_cpp';
date = '2019-10-30';
filenum = 1;

%% Load data
dirname =strcat('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data','\',date,'\',motor,'\',mt_type,'\'); %windows
% dirname =strcat('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/in_vitro_data','/',date,'/',motor,'/',mt_type,'/'); %mac

% Read in microtubule data
mt_file = dir(fullfile(dirname,'MT*.csv')); %finds appropriate file
fid=fopen(fullfile(dirname,mt_file(filenum).name)); %opens the specified file in the list and imports data
temp_mt_data = textscan(fid,'%s %s %s %s','HeaderLines',1,'Delimiter',',','EndOfLine','\n','CommentStyle','C2'); %cell with columns %1=id %2=roi_name %3=x %4=y
fclose(fid); 
[mts, interp_mts] = filter_mts(temp_mt_data, analyze_mt_num, filt_cross_mt, mt_cross_dist, pixel_size, num_pix_x, num_pix_y, zplot);
all_interp_mts = [];
for j = 1: size(interp_mts,1)
    all_interp_mts = [all_interp_mts;interp_mts{j}];
end
num_mts = size(mts,1);

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
    figure, traj_plot=gcf;
    xlabel('x (nm)'), ylabel('y (nm)'), title('Trajectories')
    
    figure, kymo_plot=gcf;
    xlabel('time (s)'), ylabel('position (nm)'), title('Kymographs')
    
    figure, offaxis_plot=gcf;
    xlabel('time (s)'), ylabel('off-axis position (nm)'), title('Off-axis position')
end
%% Initialize variables
censored = []; %track reaches end of MT
cum_run_length = [];
cum_mts = [];
ftk = 0; %number of tracks that passes all filtering

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
                uni_ind = sort(uni_ind);
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
    
%     % for testing, fit "MT" to trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     pseudo_mt = polyfit(x_tk,y_tk,1);
%     mt_ycoords = polyval(pseudo_mt,x_tk);
%     mt_coords = [x_tk, mt_ycoords];
%     [~,uni_ind] = unique(mt_coords(:,1),'stable'); %find repeated x-values in MT
%     uni_ind = sort(uni_ind);
%     mt_coords = mt_coords(uni_ind,:);
    
    % If MT is oriented the "wrong" way (i.e. motor moves right to left on screen,as all MTs are drawn left to right), reverse MT
    if x_tk(1) > x_tk(end)
        mts{mt_tk} = flipud(mts{mt_tk});
        interp_mts{mt_tk} = flipud(interp_mts{mt_tk});
    end
    
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

    interp_mt_coords = interp_mts{mt_tk};
    % project track along MT
    fit_vector = [];
    motor_vector = [0,0;diff(x_tk),diff(y_tk)];
    mt_vector = [0,0;diff(mt_coords)];
    % polynomial fit to vector (gives MT coordinates for each spot localized in track)
    [~,uni_ind] = unique(interp_mt_coords(:,1),'stable'); %find repeated x-values in MT
    %uni_ind = sort(uni_ind);
    fit_mt = polyfit(interp_mt_coords(uni_ind,1),interp_mt_coords(uni_ind,2),1);%pchip(MTs{chosen_MT}(uni_ind,1),MTs{chosen_MT}(uni_ind,2)); %%%this still "overfits" the MT
    interp_mt_for_motor = polyval(fit_mt,x_tk);
%     fit_vector(:,1) = interp1(mt_coords(:,1),mt_vector(:,1),x_tk,'pchip'); 
%     fit_vector(:,2) = interp1(mt_coords(:,1),mt_vector(:,2),x_tk,'pchip');
    fit_vector(:,1) = x_tk;
    fit_vector(:,2) = interp_mt_for_motor;
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
    
%     moveupmt = polyfit(frame_tk,position',1);
%     if moveupmt(1) < 0 %motor moves from high to low x values (i.e. MT drawn the "wrong" way)
% %         position = position.*-1;
%         
%         flip_mt = 1;
%         mts{mt_tk} = flipud(mts{mt_tk});
%         interp_mts{mt_tk} = flipud(interp_mts{mt_tk});
%         %reverse mt_coords so motor moves "up" MT
%         mt_coords = flipud(mt_coords);
%         interp_mt_coords = interp_mts{mt_tk}; %flipud(interp_mts{mt_tk});
%         
%         % project track along MT
%         fit_vector = [];
%         motor_vector = [0,0;diff(x_tk),diff(y_tk)];
%         mt_vector = [0,0;diff(mt_coords)];
%         % polynomial fit to vector (gives MT coordinates for each spot localized in track)
%         [~,uni_ind] = unique(interp_mt_coords(:,1),'stable'); %find repeated x-values in MT
%         %uni_ind = sort(uni_ind);
%         fit_mt = polyfit(interp_mt_coords(uni_ind,1),interp_mt_coords(uni_ind,2),1);%pchip(MTs{chosen_MT}(uni_ind,1),MTs{chosen_MT}(uni_ind,2)); %%%this still "overfits" the MT
%         interp_mt_for_motor = polyval(fit_mt,x_tk);
%     %     fit_vector(:,1) = interp1(mt_coords(:,1),mt_vector(:,1),x_tk,'pchip'); 
%     %     fit_vector(:,2) = interp1(mt_coords(:,1),mt_vector(:,2),x_tk,'pchip');
%         fit_vector(:,1) = x_tk;
%         fit_vector(:,2) = interp_mt_for_motor;
%         fit_norm = repmat(sqrt(sum(fit_vector.^2,2)),1,2);
%         fit_vector_norm = fit_vector./fit_norm; %normalize
%         %project
%         delp = zeros(1,numel(x_tk));
%         delp_off = zeros(1,numel(x_tk));
%         for kv = 1:numel(x_tk)
%             delp(kv) = motor_vector(kv,:)*fit_vector_norm(kv,:)'; %displacement projected onto MT axis
%             delp_off(kv) = motor_vector(kv,:)*(fit_vector_norm(kv,:)*[0,1;-1,0])'; %displacement orthogonal to MT axis
%         end
%         position = cumsum(delp); %position of track along MT vector
%         position_off = cumsum(delp_off); %position of track perpendicular to MT vector
%     else
%         flip_mt = 0; %did not flip MT
%     end

%     % Plot trajectory
%     if zplot ~= 0 %&& skip_tk~=1 %&& mod(tk,10) == 0
%         figure(traj_plot), hold on, plot(x_tk,y_tk,'.-','Color',cmap(tk,:)) %plots trajectories
%         figure(kymo_plot), hold on, plot(frame_tk.*exp_time,position,'.-','Color',cmap(tk,:)) %plots "kymographs"
%     end
    
    %% Sliding MSD analysis
    if length (frame_tk) > l_window + 4
        tmsd_res = tMSD_2D(x_tk,y_tk,frame_tk,l_window,exp_time,msd_thresh, msd_step,l_min);
        tmsd_frames = tmsd_res(:,1); %frame of particle (in movie)
        loc_alpha = tmsd_res(:,2); %local alpha-value
        proc_frames = tmsd_res(:,3); %was frame marked as processive (1) or paused/diffusive (0)
    end
    
    %%%%%%%%%%%%%% All filtering is done before this point %%%%%%%%%%%%%%%%
    ftk = ftk+1; %only advances if track is not filtered
    censored = [censored; cens_tk]; %cens_tk is 1 if track reaches MT end, 0 if not
    
    %% Calculate parameters
    inst_vel = diff(position)./diff(frame_tk.*exp_time);
    if length (frame_tk) > l_window + 4
        proc_vel = diff(position(proc_frames==1))./diff(frame_tk(proc_frames==1).*exp_time);
    end
    
    %% Store data
    traj(ftk).position = position;
%     traj(ftk).flip = flip_mt;
    traj(ftk).mt_coords = mt_coords; %will be flipped if needed for track --> does not necessarily match mts{mt_tk}
    traj(ftk).interp_mt_coords = interp_mt_coords; %will be flipped if needed for track --> does not necessarily match interp_mts{mt_tk}
    traj(ftk).offaxis_position = position_off;
    traj(ftk).x = x_tk;
    traj(ftk).y = y_tk;
    traj(ftk).frames = frame_tk;
    traj(ftk).mt = mt_tk;
    traj(ftk).mt_length = sum(sqrt(mt_vector(:,1).^2 + mt_vector(:,2).^2));
    traj(ftk).run_length = position(end);
    traj(ftk).inst_vel = inst_vel;
    if length (frame_tk) > l_window + 4
        traj(ftk).proc_frames = proc_frames;
        traj(ftk).loc_alpha = loc_alpha;
        traj(ftk).proc_vel = proc_vel;
    end
    
    cum_run_length = [cum_run_length; abs(position(end))];
    cum_censored = censored;
    cum_mts = [cum_mts; mt_tk];
    
end
% profile viewer

%% Plot
if zplot ~= 0
    for ftk = 1:size(cum_run_length,1)
        figure(traj_plot), hold on, plot(traj(ftk).x,traj(ftk).y,'.-','Color',cmap(ftk,:)) %plots trajectories
        figure(kymo_plot), hold on, plot((traj(ftk).frames).*exp_time,traj(ftk).position,'.-','Color',cmap(ftk,:)) %plots "kymographs"
        figure(offaxis_plot), hold on, plot((traj(ftk).frames).*exp_time,traj(ftk).offaxis_position,'.-','Color',cmap(ftk,:)) %plots off-axis "kymographs"
    end
    
    for mttk = 1:num_mts
        ftk_on_mt = find(cum_mts == mttk); %gives indices of cum_mts, which should match that of ftk
        if ~isempty(ftk_on_mt)
            mt_ends = [mts{mttk}(1,:);mts{mttk}(end,:)]';
%             mt_xy = interp_mts{mttk};
             mt_vector = [0,0;diff(mts{mttk})];
             mt_cum_length = cumsum(sqrt(mt_vector(:,1).^2 + mt_vector(:,2).^2));
             mt_length = sum(sqrt(mt_vector(:,1).^2 + mt_vector(:,2).^2));

            % polynomial fit to vector (gives MT coordinates for each spot localized in track)
            [~,uni_ind] = unique(interp_mts{mttk}(:,1),'stable'); %find repeated x-values in MT
            %uni_ind = sort(uni_ind);
            fit_mt = polyfit(interp_mts{mttk}(uni_ind,1),interp_mts{mttk}(uni_ind,2),1);%pchip(MTs{chosen_MT}(uni_ind,1),MTs{chosen_MT}(uni_ind,2)); %%%this still "overfits" the MT
            mt_slope = fit_mt(1);
            mt_yint = fit_mt(2);
            perp_slope = -1/mt_slope;

            figure, hold on
            for j=1:length(ftk_on_mt)
                %find closest point along MT to start of track
                perp_yint_start = -perp_slope*traj(ftk_on_mt(j)).x(1)+traj(ftk_on_mt(j)).y(1);
                x_intersect_start = (perp_yint_start-mt_yint)/(mt_slope-perp_slope);
                y_intersect_start = perp_slope*x_intersect_start + perp_yint_start;
                %find closest point along MT to end of track
                perp_yint_end = -perp_slope*traj(ftk_on_mt(j)).x(end)+traj(ftk_on_mt(j)).y(end);
                x_intersect_end = (perp_yint_end-mt_yint)/(mt_slope-perp_slope);
                y_intersect_end = perp_slope*x_intersect_end + perp_yint_end;
                
                %check is end further from MT start than start? if yes, MT
                %is the right way around add start to all, if not, add MT
                %length - start to all
                dist_to_start = sqrt((mt_ends(1,1) - x_intersect_start).^2+(mt_ends(1,2) - y_intersect_start).^2); %distance from start of MT to start of track
                dist_to_end = sqrt((mt_ends(1,1) - x_intersect_end).^2+(mt_ends(1,2) - y_intersect_end).^2); %distance from start of MT to end of track (if this is smaller, MT is the "wrong" way around
                %mt_start = sqrt(mt_ends(1,1)^2 + mt_ends(1,2)^2);
                
                if dist_to_start <= dist_to_end %MT is the right way around (motor moves "up" the MT)
                    pos_on_mt = traj(ftk_on_mt(j)).position + dist_to_start;
                else
                    %other end of MT is start
                    dist_to_start = sqrt((mt_ends(2,1) - x_intersect_start).^2+(mt_ends(2,2) - y_intersect_start).^2); %distance from start of MT to start of track
                    dist_to_end = sqrt((mt_ends(2,1) - x_intersect_end).^2+(mt_ends(2,2) - y_intersect_end).^2); %distance from start of MT to end of track (if this is smaller, MT is the "wrong" way around
                    %mt_start = sqrt(mt_ends(2,1)^2 + mt_ends(2,2)^2);
                    pos_on_mt = traj(ftk_on_mt(j)).position + dist_to_start;
                end

                plot((traj(ftk_on_mt(j)).frames).*exp_time,pos_on_mt,'.-','Color',cmap(ftk_on_mt(j),:)) %plots "kymographs"
                %plot((traj(ftk_on_mt(j)).frames).*exp_time,traj(ftk_on_mt(j)).position,'.-','Color',cmap(ftk_on_mt(j),:)) %plots "kymographs"
            end
            xlabel('time (s)'), ylabel('position along MT (nm)'), title(['Kymographs for MT number ', num2str(mttk), ' (MT length ', num2str(mt_length), 'nm)'])
            hold off
        end
    end
end


%% Plot parameters to check data
if zplot ~=0
    figure,totrl=gcf; %initialize figure
    [trl_n, trl_edges]=histcounts(cum_run_length, 'BinWidth', rl_binwidth, 'Normalization', 'pdf');
    nhist_trl=trl_n;
    xhist_trl=trl_edges+(rl_binwidth/2);
    xhist_trl(end)=[]; 
    figure(totrl), hold on 
    bar(xhist_trl,nhist_trl)
    xlabel('Total run length (nm)'), ylabel('Probability density'), title([motor,' ', mt_type,' Total run length histogram'])
    opt = statset('mlecustom');
    opt = statset(opt,'FunValCheck','off','MaxIter',1e5,'MaxFunEvals',1e5,'Display','iter','TolFun',10e-20);
    p0 = [1000];
    loL = [300];
    upL = [3000];
    [estimRL, pciRL] = mle(cum_run_length,'Distribution','exponential','Censoring',censored,'start',p0,'Options',opt,'LowerBound', loL);%'UpperBound', upL) %
    %plot fit
    yRL = exppdf(xhist_trl, estimRL);
    yRLlo = exppdf(xhist_trl, pciRL(1));
    yRLhi = exppdf(xhist_trl, pciRL(2));
    plot(xhist_trl,0.8.*yRL,'r-');
    plot(xhist_trl,0.8.*yRLlo,'r.');
    plot(xhist_trl,0.8.*yRLhi,'r.');
    hold off

    
end

%% Save data
if zsave ~= 0
    save_dirname =strcat('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data\results'); %windows
    % save_dirname =strcat('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/in_vitro_data/results'); %mac
    save_filename = ['post_particle_tracking','_',date,'_',motor,'_',mt_type,'_',num2str(filenum)];
    
    save([save_dirname,save_filename],'mts','interp_mts','traj','cum_run_length','cum_censored')
end