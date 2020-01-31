%% UU - Kapitein Lab_
% Analyze in vitro single molecule motility assays
% MKI and CC 2019-11-05

%%
% TO DO:
% MSD analysis for vel, diff coeff
% sliding MSD analysis

%% This code reads in the files generated by post_particle_tracking_1 for one condition.
% To use this file you will need the following files:
% .mat file output by post_particle_tracking_1
%
% Ensure that tMSD_2D.m and MSD_2D.m are in the same folder as this code.

clear all, close all
addpath('C:\Users\6182658\OneDrive - Universiteit Utrecht\MATLAB\GitHub Codes\in-vitro-codes\kapitein-invitro') %windows
%addpath('/Users/malinaiwanski/Documents/MATLAB/GitHub/kapitein-invitro') %mac
addpath('C:\Users\6182658\OneDrive - Universiteit Utrecht\MATLAB') %windows
set(0,'DefaultFigureWindowStyle','docked')

%% Options (make 0 to NOT perform related action, 1 to perform)
zplot = 1;
zsave = 0;
zcap = 1; %set to 1 if using capped MTs

%% Parameters
% From imaging:
pixel_size = 64.0; %pixel size of camera [nm]
exp_time = 0.1; %exposure time [s]
num_frames = 600; %number of frames in movie [frames]
num_pix_x = 512;
num_pix_y = 512;

% For analysis
l_window = 7; %number of frames to average for sliding MSD analysis
msd_thresh = 1.1; %alpha-value above which is processive, below which is paused
msd_step = 0.3; %minimum threshold for findchangepts function; minimum improvement in residual error; changepoints currently based on mean values, but can use rms, std, mean and slope
l_min = 3; %minimum distance between two changepoints - smallest duration of pause/run

% For plotting:
rl_binwidth = 200; %bin width for run length histograms
vel_binwidth = 100; %bin width for velocity histograms
time_binwidth = 0.5; %bin width for association time histograms
loca_binwidth = 0.1; %bin width for local alpha-values (MSD analysis)

%% Movies to analyze
motor = {'kif1a','kif5b'}; %{'kif1a'}; %
mt_type = {'cap','taxol_cap'}; %{'1cycle_cpp','2cycle_cpp','gdp_taxol'}; %
date = {'2019-12-09','2019-12-13'}; %{'2019-10-30'}; %

%% Initialize figures
if zplot ~= 0
    figure,totrl=gcf; %initialize figure
    figure,totuvel=gcf; %initialize figure
    figure,boundtime=gcf; %initialize figure
    figure,localalpha=gcf; %initialize figure
    figure,instvel=gcf; %initialize figure
    figure,procvel=gcf; %initialize figure
    figure,landpos = gcf;
    figure,normlandpos = gcf;
    figure, landtime = gcf;
    figure,numtracks = gcf;
    figure,numtracksmt = gcf;
    figure, timebwland = gcf;
    figure,timebwlandhist = gcf;
    figure,timebwlandbydisthist = gcf;
    figure, landingdeltime = gcf;
    figure, landingdeltimeviolin = gcf;
%     figure, landrelativemotor = gcf;
%     figure, landrelativemotorviolin = gcf;
%     figure, ldrelmot = gcf;
    figure, mtlengthfrommotor = gcf;
    figure, totmtlength = gcf;
    figure, landingdistancetomotor = gcf;
%     figure, normlandingdistancetomotor = gcf;
%     figure, normmtlandingdistancetomotor = gcf;
    if zcap == 1
        figure,segmentvel=gcf; %initialize figure 
        figure,segmenttypevel=gcf; %initialize figure 
        figure,trackstartsegment=gcf; %initialize figure
        figure, trackstartsegmenttype=gcf;
    end
end
%% Initialize variables

%% Load data
dirname =strcat('C:\Users\6182658\OneDrive - Universiteit Utrecht\in_vitro_data\results'); %windows
%dirname =strcat('/Users/malinaiwanski/OneDrive - Universiteit Utrecht/in_vitro_data/results'); %mac
filename_start = 'post_particle_tracking';

num_cat = size(motor,2)*size(mt_type,2);
catk = 0;

for mk = 1:size(motor,2)
    for mtk = 1:size(mt_type,2)
        catk = catk+1;
        
        datcat(catk).name = strcat(string(motor(mk)),'_',string(mt_type(mtk)));
        %traj_files = dir(fullfile(dirname,[filename_start,'*',motor(mk),'_',mt_type(mtk),'*.mat']));
        traj_files = dir(fullfile(dirname,strcat('*',datcat(catk).name,'*.mat')));
        
        disp('------------Analyzing------------')
        disp(strcat('category: ',num2str(catk)))
        disp(strcat('motor: ',string(motor(mk))))
        disp(strcat('mt_type: ',string(mt_type(mtk))))
        disp(strcat('number of movies: ',num2str(numel(traj_files))))
        disp('---------------------------------')
        
        
        %% initialize variables
        datcat(catk).traj = {};
        datcat(catk).mt_id = {};
        datcat(catk).mts = {};
        datcat(catk).interp_mts = {};
        datcat(catk).mt_lengths = {};
        if zcap == 1
            datcat(catk).boundaries_on_mt = {};
            datcat(catk).segment_lengths = {};
        end
        datcat(catk).inst_vel = {};
        datcat(catk).proc_vel = {};
        %datcat(catk).loc_alpha = {};
        %datcat(catk).pos_on_mt = {};
        datcat(catk).track_start_times = {};
        
        datcat(catk).cum_run_length = [];
        datcat(catk).cum_association_time = [];
        datcat(catk).cum_censored = [];
        datcat(catk).cum_mean_vel = [];
        datcat(catk).cum_inst_vel = [];
        datcat(catk).cum_proc_vel = [];
        datcat(catk).cum_loc_alpha = [];
        datcat(catk).cum_norm_landing_pos = [];
        datcat(catk).cum_landing_dist_to_mt_end = [];
        datcat(catk).all_landing_dist = [];
%         datcat(catk).all_norm_landing_dist = [];
%         datcat(catk).all_mt_landing_dist = [];
        datcat(catk).all_dist_to_plus = [];
        datcat(catk).all_dist_to_minus = [];
        datcat(catk).all_time_diff_bw_land = [];
        datcat(catk).time_diff_bw_land = {};
        if zcap == 1
            datcat(catk).cum_plus_cap_vel = [];
            datcat(catk).cum_plus_gdp_vel = [];
            datcat(catk).cum_seed_vel = [];
            datcat(catk).cum_minus_cap_vel = [];
            datcat(catk).cum_minus_gdp_vel = [];
            datcat(catk).cum_track_start_segment = [];
        end
        
        cum_time_bw_landing = [];
        cum_time_bw_landing_by_length = [];
        cum_segment_lengths = zeros(6,1); %zeros(numel(traj_files),1);
        
        %% read in .mat file from each movie
        for movk = 1:numel(traj_files)
            %read in file
            datmovk = load(fullfile(dirname,traj_files(movk).name));
            
            disp('------------Reading in data from file:------------')
            disp(traj_files(movk).name)
            
            if numel(datmovk.cum_run_length) ~= numel(datmovk.traj)
                disp('ERROR: Data mismatch - try re-running post_particle_tracking_1')
            end
            
            %save data
            datcat(catk).traj{movk} = datmovk.traj;
            datcat(catk).mt_id{movk} = datmovk.traj.mt;
            datcat(catk).mts{movk} = datmovk.mts;
            datcat(catk).interp_mts{movk} = datmovk.interp_mts;
            datcat(catk).mt_lengths{movk} = datmovk.mt_lengths;
            if zcap == 1
                datcat(catk).boundaries_on_mt{movk} = datmovk.boundaries_on_mt;
                datcat(catk).segment_lengths{movk} = datmovk.segment_lengths;
                for i = 1:size(datcat(catk).segment_lengths{movk},2) %each MT
                    nzsegment_lengths = nonzeros(datcat(catk).segment_lengths{1,movk}{1,i});
                    for j = 1:size(nzsegment_lengths,1) %each segment in MT
                        cum_segment_lengths(j) = cum_segment_lengths(j) + nzsegment_lengths(j,1);
                    end
                end
            end
            
            temp_a = exist('datmovk.traj.proc_vel');
            if temp_a ~=0
                datcat(catk).inst_vel{movk} = datmovk.traj.inst_vel;
                datcat(catk).proc_vel{movk} = datmovk.traj.proc_vel;
                % datcat(catk).loc_alpha{movk} = datmovk.traj.loc_alpha;
            end
            
            %datcat(catk).pos_on_mt{movk} = datmovk.traj.position_on_mt;
            datcat(catk).track_start_times{movk} = datmovk.track_start_times;
            
            datcat(catk).cum_run_length = [datcat(catk).cum_run_length; datmovk.cum_run_length];
            datcat(catk).cum_censored = [datcat(catk).cum_censored; datmovk.cum_censored];
            datcat(catk).cum_association_time = [datcat(catk).cum_association_time; datmovk.cum_association_time];
            datcat(catk).cum_mean_vel = [datcat(catk).cum_mean_vel; datmovk.cum_mean_vel];
            datcat(catk).cum_inst_vel = [datcat(catk).cum_inst_vel, datmovk.cum_inst_vel]; %traj.inst_vel];
            datcat(catk).cum_proc_vel = [datcat(catk).cum_proc_vel, datmovk.cum_proc_vel]; %traj.proc_vel];
            datcat(catk).cum_loc_alpha = [datcat(catk).cum_loc_alpha;  datmovk.cum_loc_alpha];
            datcat(catk).cum_norm_landing_pos = [datcat(catk).cum_norm_landing_pos; datmovk.cum_norm_landing_pos];
            datcat(catk).cum_landing_dist_to_mt_end = [datcat(catk).cum_landing_dist_to_mt_end; datmovk.cum_landing_dist_to_mt_end];
            datcat(catk).all_landing_dist = [datcat(catk).all_landing_dist; datmovk.all_landing_dist];
%             datcat(catk).all_norm_landing_dist = [datcat(catk).all_norm_landing_dist; datmovk.all_norm_landing_dist];
%             datcat(catk).all_mt_landing_dist = [datcat(catk).all_mt_landing_dist; datmovk.all_mt_landing_dist];
            datcat(catk).all_dist_to_minus = [datcat(catk).all_dist_to_minus; datmovk.all_dist_to_minus];
            datcat(catk).all_dist_to_plus = [datcat(catk).all_dist_to_plus; datmovk.all_dist_to_plus];
            datcat(catk).all_time_diff_bw_land = [datcat(catk).all_time_diff_bw_land; datmovk.all_time_diff_bw_land];
            datcat(catk).time_diff_bw_land{movk} = datmovk.time_diff_bw_land;
            if zcap == 1
                datcat(catk).cum_plus_cap_vel = [datcat(catk).cum_plus_cap_vel; datmovk.cum_plus_cap_vel'];
                datcat(catk).cum_plus_gdp_vel = [datcat(catk).cum_plus_gdp_vel; datmovk.cum_plus_gdp_vel'];
                datcat(catk).cum_seed_vel = [datcat(catk).cum_seed_vel; datmovk.cum_seed_vel'];
                datcat(catk).cum_minus_cap_vel = [datcat(catk).cum_minus_cap_vel; datmovk.cum_minus_cap_vel'];
                datcat(catk).cum_minus_gdp_vel = [datcat(catk).cum_minus_gdp_vel; datmovk.cum_minus_gdp_vel'];
                datcat(catk).cum_track_start_segment = [datcat(catk).cum_track_start_segment; datmovk.cum_track_start_segment'];
            end
            
        end
        
        if isempty(datcat(catk).cum_run_length) == 1
            continue
        end
        %% Analyze
        
        
        %% Plot
        
        %run length   
        [trl_n, trl_edges]=histcounts(datcat(catk).cum_run_length, 'BinWidth', rl_binwidth, 'Normalization', 'pdf');
        nhist_trl=trl_n;
        xhist_trl=trl_edges+(rl_binwidth/2);
        xhist_trl(end)=[]; 
        figure(totrl)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_trl,nhist_trl)
        hold on 
        xlabel('Total run length (nm)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Total run length histogram'])
        opt = statset('mlecustom');
        opt = statset(opt,'FunValCheck','off','MaxIter',1e5,'MaxFunEvals',1e5,'Display','iter','TolFun',10e-20);
        p0 = [2000];
        loL = [300];
        upL = [5000];
        [estimRL, pciRL] = mle(datcat(catk).cum_run_length(datcat(catk).cum_run_length>700),'Distribution','exponential','start',p0,'Options',opt,'LowerBound', loL);%'UpperBound', upL) %'Censoring',datcat(catk).cum_censored(datcat(catk).cum_run_length>700)
        %plot fit
        yRL = exppdf(xhist_trl, estimRL);
        yRLlo = exppdf(xhist_trl, pciRL(1));
        yRLhi = exppdf(xhist_trl, pciRL(2));
        plot(xhist_trl,1.0.*yRL,'r-');
        plot(xhist_trl,1.0.*yRLlo,'r.');
        plot(xhist_trl,1.0.*yRLhi,'r.');
        hold off
        
        %mean velocity
        [uvel_n, uvel_edges]=histcounts(datcat(catk).cum_mean_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
        nhist_uvel=uvel_n;
        xhist_uvel=uvel_edges+(vel_binwidth/2);
        xhist_uvel(end)=[]; 
        figure(totuvel)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_uvel,nhist_uvel)
        hold on 
        xlabel('Mean velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Mean velocity histogram'])
        opt = statset('mlecustom');
        opt = statset(opt,'FunValCheck','off','MaxIter',1e5,'MaxFunEvals',1e5,'Display','iter','TolFun',10e-20);
        p0 = [200, 200];
        loL = [0, 0];
        upL = [1000, 1000];
        [estimuV, pciuV] = mle(datcat(catk).cum_mean_vel,'Distribution','normal','start',p0,'Options',opt,'LowerBound', loL);%'UpperBound', upL) %
        %plot fit
        yuV = normpdf(xhist_uvel, estimuV(1), estimuV(2));
        yuVlo = normpdf(xhist_uvel, pciuV(1,1), pciuV(1,2));
        yuVhi = normpdf(xhist_uvel, pciuV(2,1), pciuV(2,2));
        plot(xhist_uvel,1.0.*yuV,'r-');
        plot(xhist_uvel,1.0.*yuVlo,'r.');
        plot(xhist_uvel,1.0.*yuVhi,'r.');
        hold off

        %association time
        [bt_n, bt_edges]=histcounts(datcat(catk).cum_association_time, 'BinWidth', time_binwidth, 'Normalization', 'pdf');
        nhist_bt=bt_n;
        xhist_bt=bt_edges+(time_binwidth/2);
        xhist_bt(end)=[]; 
        figure(boundtime)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_bt,nhist_bt)
        hold on 
        xlabel('Association time (s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Track duration histogram'])
        opt = statset('mlecustom');
        opt = statset(opt,'FunValCheck','off','MaxIter',1e5,'MaxFunEvals',1e5,'Display','iter','TolFun',10e-20);
        p0 = [5];
        loL = [1];
        upL = [50];
        [estimbt, pcibt] = mle(datcat(catk).cum_association_time,'Distribution','exponential','start',p0,'Options',opt,'LowerBound', loL);%'UpperBound', upL) %'Censoring',datcat(catk).cum_censored,
        %plot fit
        ybt = exppdf(xhist_bt, estimbt);
        ybtlo = exppdf(xhist_bt, pcibt(1));
        ybthi = exppdf(xhist_bt, pcibt(2));
        plot(xhist_bt,1.0.*ybt,'r-');
        plot(xhist_bt,1.0.*ybtlo,'r.');
        plot(xhist_bt,1.0.*ybthi,'r.');
        hold off
        
        %local alpha-value
        [loca_n, loca_edges]=histcounts(datcat(catk).cum_loc_alpha, 'BinWidth', loca_binwidth, 'Normalization', 'pdf');
        nhist_loca=loca_n;
        xhist_loca=loca_edges+(loca_binwidth/2);
        xhist_loca(end)=[]; 
        opt = statset('MaxIter',1e5,'Display','iter','TolFun',1e-20);
        sigmastart = cat(3, 0.01, 0.4);
        p0 = struct('mu', [0; 2], 'Sigma', sigmastart, 'ComponentProportion', [0.5; 0.5]);
        gmm_loc_alpha = fitgmdist(datcat(catk).cum_loc_alpha, 2, 'start', p0, 'options', opt); %
        y_loca = pdf(gmm_loc_alpha,xhist_loca');
        figure(localalpha)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_loca,nhist_loca)
        hold on 
        plot(xhist_loca,y_loca,'-')
        xlabel('Local alpha-value'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Local alpha-value histogram'])
        hold off
        
        %instantaneous velocity
        [instvel_n, instvel_edges]=histcounts(datcat(catk).cum_inst_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
        nhist_instvel=instvel_n;
        xhist_instvel=instvel_edges+(vel_binwidth/2);
        xhist_instvel(end)=[]; 
        opt = statset('MaxIter',1e5,'Display','iter','TolFun',1e-20);
        sigmastart = cat(3, 0.1, 0.1);
        p0 = struct('mu', [0; 1500], 'Sigma', sigmastart, 'ComponentProportion', [0.2; 0.8]);
        gmm_inst_vel = fitgmdist(datcat(catk).cum_inst_vel', 2, 'options', opt); %'start', p0, 
        y_instvel = pdf(gmm_inst_vel,xhist_instvel');
        figure(instvel)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_instvel,nhist_instvel)
        hold on 
        plot(xhist_instvel,y_instvel,'-')
        xlabel('Instantaneous velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Instantaneous velocity histogram'])
        hold off
        
        %instantaneous processive velocity
        [procvel_n, procvel_edges]=histcounts(datcat(catk).cum_proc_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
        nhist_procvel=procvel_n;
        xhist_procvel=procvel_edges+(vel_binwidth/2);
        xhist_procvel(end)=[]; 
        opt = statset('mlecustom');
        opt = statset(opt,'FunValCheck','off','MaxIter',1e5,'MaxFunEvals',1e5,'Display','iter','TolFun',10e-20);
        p0 = [200, 200];
        loL = [0, 0];
        upL = [1000, 1000];
        [estimprocV, pciprocV] = mle(datcat(catk).cum_proc_vel,'Distribution','normal','start',p0,'Options',opt,'LowerBound', loL);%'UpperBound', upL) %
        figure(procvel)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_procvel,nhist_procvel)
        hold on 
        xlabel('Processive velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Processive velocity histogram'])
        yprocV = normpdf(xhist_procvel, estimprocV(1), estimprocV(2));
        yprocVlo = normpdf(xhist_procvel, pciprocV(1,1), pciprocV(1,2));
        yprocVhi = normpdf(xhist_procvel, pciprocV(2,1), pciprocV(2,2));
        plot(xhist_procvel,1.0.*yprocV,'r-');
        plot(xhist_procvel,1.0.*yprocVlo,'r.');
        plot(xhist_procvel,1.0.*yprocVhi,'r.');
        hold off
        
        % velocity on different segments
        if zcap ==1
            %each segment
            [seg1vel_n, seg1vel_edges]=histcounts(datcat(catk).cum_plus_cap_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_seg1vel=seg1vel_n;
            xhist_seg1vel=seg1vel_edges+(vel_binwidth/2);
            xhist_seg1vel(end)=[]; 
            [seg2vel_n, seg2vel_edges]=histcounts(datcat(catk).cum_plus_gdp_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_seg2vel=seg2vel_n;
            xhist_seg2vel=seg2vel_edges+(vel_binwidth/2);
            xhist_seg2vel(end)=[]; 
            [seg3vel_n, seg3vel_edges]=histcounts(datcat(catk).cum_seed_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_seg3vel=seg3vel_n;
            xhist_seg3vel=seg3vel_edges+(vel_binwidth/2);
            xhist_seg3vel(end)=[]; 
            [seg4vel_n, seg4vel_edges]=histcounts(datcat(catk).cum_minus_gdp_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_seg4vel=seg4vel_n;
            xhist_seg4vel=seg4vel_edges+(vel_binwidth/2);
            xhist_seg4vel(end)=[]; 
            [seg5vel_n, seg5vel_edges]=histcounts(datcat(catk).cum_minus_cap_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_seg5vel=seg5vel_n;
            xhist_seg5vel=seg5vel_edges+(vel_binwidth/2);
            xhist_seg5vel(end)=[]; 
            figure(segmentvel) 
            subplot(size(motor,2),size(mt_type,2),catk)
            bar(xhist_seg1vel,nhist_seg1vel,'m','FaceAlpha',0.6)
            hold on 
            bar(xhist_seg2vel,nhist_seg2vel,'g','FaceAlpha',0.6)
            bar(xhist_seg3vel,nhist_seg3vel,'k','FaceAlpha',0.6)
            bar(xhist_seg4vel,nhist_seg4vel,'b','FaceAlpha',0.6)
            bar(xhist_seg5vel,nhist_seg5vel,'r','FaceAlpha',0.6)      
            xlabel('Instantaneous velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Instantaneous velocity by segment histogram'])
            hold off
            
            %segment type
            [capvel_n, capvel_edges]=histcounts([datcat(catk).cum_plus_cap_vel;datcat(catk).cum_minus_cap_vel], 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_capvel=capvel_n;
            xhist_capvel=capvel_edges+(vel_binwidth/2);
            xhist_capvel(end)=[];
            [gdpvel_n, gdpvel_edges]=histcounts([datcat(catk).cum_plus_gdp_vel;datcat(catk).cum_minus_gdp_vel], 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_gdpvel=gdpvel_n;
            xhist_gdpvel=gdpvel_edges+(vel_binwidth/2);
            xhist_gdpvel(end)=[];
            [seedvel_n, seedvel_edges]=histcounts(datcat(catk).cum_seed_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
            nhist_seedvel=seedvel_n;
            xhist_seedvel=seedvel_edges+(vel_binwidth/2);
            xhist_seedvel(end)=[];
            figure(segmenttypevel)
            subplot(size(motor,2),size(mt_type,2),catk)
            bar(xhist_capvel,nhist_capvel,'m','FaceAlpha',0.6)
            hold on 
            bar(xhist_gdpvel,nhist_gdpvel,'g','FaceAlpha',0.6)
            bar(xhist_seedvel,nhist_seedvel,'k','FaceAlpha',0.6)     
            xlabel('Instantaneous velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Instantaneous velocity by segment type histogram'])
            hold off
        end
        
        % landing position along MT - distance from plus-end
        [landpos_n, landpos_edges]=histcounts(datcat(catk).cum_landing_dist_to_mt_end, 'BinWidth', 500, 'Normalization', 'count');
        nhist_landpos=landpos_n;
        xhist_landpos=landpos_edges+(500/2);
        xhist_landpos(end)=[]; 
        figure(landpos)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_landpos,nhist_landpos)
        hold on 
        xlabel('Landing distance to MT plus-end (nm)'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Landing distance from MT plus-end'])
        hold off
        
        % landing position along MT - normalized
        [normlandpos_n, normlandpos_edges]=histcounts(datcat(catk).cum_norm_landing_pos, 'BinWidth', 0.1, 'Normalization', 'count');
        nhist_normlandpos=normlandpos_n;
        xhist_normlandpos=normlandpos_edges+(0.1/2);
        xhist_normlandpos(end)=[]; 
        figure(normlandpos)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_normlandpos,nhist_normlandpos)
        hold on 
        xlabel('Landing position along MT (fraction of 1)'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Normalized Landing Position'])
        hold off
        
        % landing time on given MT
        figure(landtime)
        subplot(size(motor,2),size(mt_type,2),catk)
        num_tracks_mt = [];
        mt_lengths = [];
        num_tracks_per_um = [];
        num_mov = size(datcat(catk).track_start_times,2);
        mtnum = 0;
        for movj = 1:num_mov
            num_mt_mov = size(datcat(catk).track_start_times{1,movj},1);
            for mtj = 1:num_mt_mov
                if isempty(datcat(catk).mts{movj}{mtj})
                    continue
                else
                    mtnum = mtnum+1;
                    mt_length_temp = (datcat(catk).mt_lengths{movj}{mtj})/1000;
    %                 if ~isempty(mt_length_temp)
                        mt_lengths(mtnum) = mt_length_temp;
    %                 end
                    num_tracks_mt(mtnum) = size(cell2mat(datcat(catk).track_start_times{1,movj}(mtj,1)),1);
                    num_tracks_per_um_temp = num_tracks_mt(mtnum)/mt_length_temp;
                    num_tracks_per_um = [num_tracks_per_um; num_tracks_per_um_temp];
                    cmap = colormap(parula(num_mov));
                    scatter(repmat(mtnum,[num_tracks_mt(mtnum),1]), cell2mat(datcat(catk).track_start_times{1,movj}(mtj,1)).*mt_lengths(mtnum),25,cmap(movj,:),'filled')
                    hold on 
                end
            end
        end
        xlabel('MT number'), ylabel('Track start time (s*\mum)'), title([motor{mk},' ', mt_type{mtk},' Normalized Landing times on each MT'])
%         xlim([0 mtnum])
        hold off
        % number of tracks on a MT per um
        [numtkpum_n, numtkpum_edges]=histcounts(num_tracks_per_um, 'BinWidth', 0.5, 'Normalization', 'count');
        nhist_numtkpum=numtkpum_n;
        xhist_numtkpum=numtkpum_edges+(0.5/2);
        xhist_numtkpum(end)=[]; 
        figure(numtracks)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_numtkpum,nhist_numtkpum)
        hold on 
        xlabel('Number of track starts per \mum'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Normalized Number of Tracks on given MT'])
        hold off
        % number of tracks on a MT
        [numtkpmt_n, numtkpmt_edges]=histcounts(num_tracks_mt, 'BinWidth', 2, 'Normalization', 'probability');
        nhist_numtkpmt=numtkpmt_n;
        xhist_numtkpmt=numtkpmt_edges+(2/2);
        xhist_numtkpmt(end)=[]; 
        figure(numtracksmt)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_numtkpmt,nhist_numtkpmt)
        hold on 
        xlabel('Number of track starts per MT'), ylabel('Probability'), title([motor{mk},' ', mt_type{mtk},' Number of Tracks on given MT'])
        hold off
        
        %time between landing events on a given MT
        figure(timebwland)
        subplot(size(motor,2),size(mt_type,2),catk)
        num_mov = size(datcat(catk).track_start_times,2);
        mtnum = 0;
        for movj = 1:num_mov
            num_mt_mov = size(datcat(catk).track_start_times{1,movj},1);
            for mtj = 1:num_mt_mov
                if isempty(datcat(catk).mts{movj}{mtj})
                    continue
                else
                    mtnum = mtnum+1;
                    num_tracks_mt = size(cell2mat(datcat(catk).track_start_times{1,movj}(mtj,1)),1);
                    mt_lengths(mtnum) = datcat(catk).mt_lengths{movj}{mtj}/1000;
                    time_bw_landing = diff(cell2mat(datcat(catk).track_start_times{1,movj}(mtj,1)));
                    cum_time_bw_landing = [cum_time_bw_landing; time_bw_landing];
                    cum_time_bw_landing_by_length = [cum_time_bw_landing_by_length; time_bw_landing*mt_lengths(mtnum)];
                    if numel(time_bw_landing) > 2
                        [idx,C,sumd] = kmeans(time_bw_landing, 2);
    %                     if C(1) < C(2)
    %                         plot(repmat(1,[size(time_bw_landing(idx==1,1),1),1]),time_bw_landing(idx==1,1),'r.','MarkerSize',12)
    %                         hold on
    %                         plot(repmat(1,[size(time_bw_landing(idx==2,1),1),1]),time_bw_landing(idx==2,1),'b.','MarkerSize',12)
    %                         %plot(C(:,1),'kx','MarkerSize',15,'LineWidth',3) 
    %                         %legend('Cluster 1','Cluster 2','Centroids', 'Location','NW')
    %                     else
    %                         plot(repmat(1,[size(time_bw_landing(idx==1,1),1),1]),time_bw_landing(idx==1,1),'b.','MarkerSize',12)
    %                         hold on
    %                         plot(repmat(1,[size(time_bw_landing(idx==2,1),1),1]),time_bw_landing(idx==2,1),'r.','MarkerSize',12)
    %                    
    %                     end
                        if C(1) < C(2)
                            plot(repmat(mtnum,[size(time_bw_landing(idx==1,1),1),1]),time_bw_landing(idx==1,1).*mt_lengths(mtnum),'r.','MarkerSize',12)
                            hold on
                            plot(repmat(mtnum,[size(time_bw_landing(idx==2,1),1),1]),time_bw_landing(idx==2,1).*mt_lengths(mtnum),'b.','MarkerSize',12)
                            %plot(C(:,1),'kx','MarkerSize',15,'LineWidth',3) 
                            %legend('Cluster 1','Cluster 2','Centroids', 'Location','NW')
                        else
                            plot(repmat(mtnum,[size(time_bw_landing(idx==1,1),1),1]),time_bw_landing(idx==1,1).*mt_lengths(mtnum),'b.','MarkerSize',12)
                            hold on
                            plot(repmat(mtnum,[size(time_bw_landing(idx==2,1),1),1]),time_bw_landing(idx==2,1).*mt_lengths(mtnum),'r.','MarkerSize',12)

                        end
                        xlabel('MT number'), ylabel('Time between track start events (s*\mum)'), title([motor{mk},' ', mt_type{mtk},' Normalized Time between start of tracks'])
                        xlim([0 mtnum])
                        %hold off
                    end
                end
            end 
        end
        [deltland_n, deltland_edges]=histcounts(cum_time_bw_landing, 'BinWidth', 0.5, 'Normalization', 'count');
        nhist_deltland=deltland_n;
        xhist_deltland=deltland_edges+(0.5/2);
        xhist_deltland(end)=[]; 
        figure(timebwlandhist)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_deltland,nhist_deltland)
        hold on 
        xlabel('Time between landing events (s)'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Time between landing events'])
        hold off
        
        [deltlanddist_n, deltlanddist_edges]=histcounts(cum_time_bw_landing_by_length, 'BinWidth', 5, 'Normalization', 'pdf');
        nhist_deltlanddist=deltlanddist_n;
        xhist_deltlanddist=deltlanddist_edges+(5/2);
        xhist_deltlanddist(end)=[]; 
        figure(timebwlandbydisthist)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_deltlanddist,nhist_deltlanddist)
        hold on 
        xlabel('Time between landing events (s*\mum)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Normalized Time between landing events'])
        hold off
         
%         figure(landrelativemotorviolin)
%         subplot(size(motor,2),size(mt_type,2),catk)
%         violinplot(datcat(catk).tot_nhist_landdist,datcat(catk).tot_xhist_landdist)
%         hold on 
%         xlabel('Distance from motor (\mum)'), ylabel('Landing rate (/\mum /s)'), title([motor{mk},' ', mt_type{mtk},' Landing analysis'])
%         hold off

        %MT length on either side of motor
        datcat(catk).all_dist_to_minus = -1.*datcat(catk).all_dist_to_minus;
        [mtlength_n, mtlength_edges]=histcounts([datcat(catk).all_dist_to_plus;datcat(catk).all_dist_to_minus], 'BinWidth', 500, 'Normalization', 'probability');
        nhist_mtlength=mtlength_n;
        landing_distances_edges = mtlength_edges;
        xhist_mtlength=mtlength_edges+(500/2);
        xhist_mtlength(end)=[]; 
        figure(mtlengthfrommotor)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_mtlength,nhist_mtlength)
        hold on 
        xlabel('Length of MT away from motor (nm)'), ylabel('Fraction'), title([motor{mk},' ', mt_type{mtk},' Length of MT on either side of motor'])
        hold off
        
        % total MT lengths
        [totmtlength_n, totmtlength_edges]=histcounts(mt_lengths, 'BinWidth', 1, 'Normalization', 'probability');
        nhist_totmtlength=totmtlength_n;
        xhist_totmtlength=totmtlength_edges+(1/2);
        xhist_totmtlength(end)=[]; 
        figure(totmtlength)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_totmtlength,nhist_totmtlength)
        hold on 
        xlabel('Length of MT (\mum)'), ylabel('Fraction'), title([motor{mk},' ', mt_type{mtk},' Length of MTs'])
        hold off
        
        %landing distance to existing track
        %[landing_distances_n, landing_distances_edges]=histcounts(datcat(catk).all_landing_dist, 'BinWidth', 500, 'Normalization', 'count');
        [landing_distances_n, landing_distances_edges]=histcounts(datcat(catk).all_landing_dist, landing_distances_edges, 'Normalization', 'count');
        nhist_landing_distances=landing_distances_n;
        xhist_landing_distances=landing_distances_edges+(500/2);
        xhist_landing_distances(end)=[]; 
        figure(landingdistancetomotor)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_landing_distances,nhist_landing_distances)
        hold on 
        xlabel('Landing distance from existing track (nm)'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Landing distance from existing track']), hold on
        xlim([-3500, 3500])
        hold off
        
        %time difference between landings
        [landing_delt_n, landing_delt_edges]=histcounts(datcat(catk).all_time_diff_bw_land, 'BinWidth', 1, 'Normalization', 'count');
        nhist_landing_delt=landing_delt_n;
        xhist_landing_delt=landing_delt_edges+(1/2);
        xhist_landing_delt(end)=[]; 
        figure(landingdeltime)
        subplot(size(motor,2),size(mt_type,2),catk)
        bar(xhist_landing_delt,nhist_landing_delt)
        hold on 
        xlabel('Time between landings (s)'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Time between landings within 6\mum'])
        xlim([0 60])
        hold off
        
        delta_t = [];
        mt_cat = [];
        mt_plot_ind = 0;
        for i = 1:size(datcat(catk).time_diff_bw_land,2)
            for j = 1: size(datcat(catk).time_diff_bw_land{1,i},2)
                mt_plot_ind = mt_plot_ind +1;
                delta_t = [delta_t; datcat(catk).time_diff_bw_land{1,i}{1,j}];
                mt_cat = [mt_cat; repmat(mt_plot_ind,size(datcat(catk).time_diff_bw_land{1,i}{1,j},1),1)];
            end
        end
%         figure(landingdeltimeviolin)
%         subplot(size(motor,2),size(mt_type,2),catk)
%         plot(mt_cat,delta_t,'.','MarkerSize',12)
% %         violinplot(delta_t,mt_cat)
%         hold on 
%         xlabel('MT'), ylabel('Time between landings (s)'), title([motor{mk},' ', mt_type{mtk},' Time between landing in first 6\mum of MT'])
%         hold off
        
%         %normalized landing distance to existing track
%         [norm_landing_distances_n, norm_landing_distances_edges]=histcounts(datcat(catk).all_norm_landing_dist, 'BinWidth', 0.1, 'Normalization', 'count');
%         nhist_norm_landing_distances=norm_landing_distances_n;
%         xhist_norm_landing_distances=norm_landing_distances_edges+(0.1/2);
%         xhist_norm_landing_distances(end)=[]; 
%         figure(normlandingdistancetomotor)
%         subplot(size(motor,2),size(mt_type,2),catk)
%         bar(xhist_norm_landing_distances,nhist_norm_landing_distances)
%         hold on 
%         xlabel('Normalized landing distance from existing track'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Normalized landing distance from existing track (MT side with motor)'])
%         hold off
%         
%         %normalized landing distance to existing track
%         [norm_mt_landing_distances_n, norm_mt_landing_distances_edges]=histcounts(datcat(catk).all_mt_landing_dist, 'BinWidth', 0.1, 'Normalization', 'count');
%         nhist_norm_mt_landing_distances=norm_mt_landing_distances_n;
%         xhist_norm_mt_landing_distances=norm_mt_landing_distances_edges+(0.1/2);
%         xhist_norm_mt_landing_distances(end)=[]; 
%         figure(normmtlandingdistancetomotor)
%         subplot(size(motor,2),size(mt_type,2),catk)
%         bar(xhist_norm_mt_landing_distances,nhist_norm_mt_landing_distances)
%         hold on 
%         xlabel('Normalized landing distance from existing track'), ylabel('Count'), title([motor{mk},' ', mt_type{mtk},' Normalized landing distance from existing track (full MT)'])
%         hold off

        %track start segment
        if zcap == 1
            [tkstartseg_n, tkstartseg_edges]=histcounts(datcat(catk).cum_track_start_segment, 'BinEdges', [0.5,1.5,2.5,3.5,4.5,5.5], 'Normalization', 'count');
            nhist_tkstartseg=tkstartseg_n;
            xhist_tkstartseg= [1,2,3,4,5];
            %xhist_tkstartseg(end)=[]; 
            for j = 1:size(xhist_tkstartseg,2)
                nhist_tkstartseg(j) = nhist_tkstartseg(j)/(cum_segment_lengths(j)/1000);
            end
            figure(trackstartsegment) 
            subplot(size(motor,2),size(mt_type,2),catk)
            bar(xhist_tkstartseg,nhist_tkstartseg)
            hold on 
            set(gca,'XTick',[1,2,3,4,5], 'xticklabel',{'Plus GMP-CPP cap','Plus GDP lattice','GMP-CPP seed','Minus GDP lattice','Minus GMP-CPP cap'});
            xlabel('MT Segment'), ylabel('Number of tracks starting in segment per unit distance (/\mum)'), title([motor{mk},' ', mt_type{mtk},' Track Start Segment'])

            nhist_tk_starttype = zeros(3,1);
            xhist_tkstarttype = [1,2,3];
            figure(trackstartsegmenttype)
            nhist_tkstarttype(1) = nhist_tkstartseg(1)+nhist_tkstartseg(5);
            nhist_tkstarttype(2) = nhist_tkstartseg(2)+nhist_tkstartseg(4);
            nhist_tkstarttype(3) = nhist_tkstartseg(3);
            subplot(size(motor,2),size(mt_type,2),catk)
            bar(nhist_tkstarttype)
            hold on 
            set(gca,'XTick',[1,2,3], 'xticklabel',{'GMP-CPP caps','GDP lattice','GMP-CPP seed'});
            xlabel('MT Segment Type'), ylabel('Number of tracks starting in segment per unit distance (/\mum)'), title([motor{mk},' ', mt_type{mtk},' Track Start Segment'])
        end 
%         failed = lifetimes(find(censored==0));
%         nfailed = size(failed);
%         survived = lifetimes(find(censored~=0));
%         nsurvived = size(survived);
%         frac_survive = nsurvived/(nfailed+nsurvived); %proportion of censored data
% 
%         figure, hold on %plot observed cluster persistence times
%         plot([zeros(nsurvived),survived]', repmat(1:nsurvived,2,1),'Color','b','LineStyle','-') %censored times in blue
%         plot([zeros(nfailed),failed]', repmat(nsurvived+1:nsurvived+nfailed,2,1),'Color','r','LineStyle','-') %others in red
%         xlabel('Survival time'); ylabel('Observation number')
%         hold off

        
    end
end

%% Analyze




%% Save data
if zsave ~= 0
    
end