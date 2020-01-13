%% UU - Kapitein Lab_
% Analyze in vitro single molecule motility assays
% MKI and CC 2019-11-05

%%
% TO DO:
% MSD analysis for vel, diff coeff
% sliding MSD analysis
% normalize landing time with MT length

%% This code reads in the files generated by post_particle_tracking_1 for one condition.
% To use this file you will need the following files:
% .mat file output by post_particle_tracking_1
%
% Ensure that tMSD_2D.m and MSD_2D.m are in the same folder as this code.

clear all, close all
addpath('C:\Users\6182658\OneDrive - Universiteit Utrecht\MATLAB\GitHub Codes\in-vitro-codes\kapitein-invitro') %windows
%addpath('/Users/malinaiwanski/Documents/MATLAB/GitHub/kapitein-invitro') %mac
set(0,'DefaultFigureWindowStyle','docked')

%% Options (make 0 to NOT perform related action, 1 to perform)
zplot = 1;
zsave = 0;
zcap = 0; %set to 1 if using capped MTs

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
mt_type = {'1cycle_cpp','2cycle_cpp','gdp_taxol'}; %{'cap'}; %
date = {'2019-10-30'}; %{'2019-12-09'}; %

%% Initialize figures
if zplot ~= 0
    
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
        cum_segment_lengths = zeros(5,1);
        
        %% read in .mat file from each movie
        for movk = 1:numel(traj_files)
            %read in file
            datmovk = load(fullfile(dirname,traj_files(movk).name));

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
                for i = 1:size(datcat(catk).segment_lengths{movk},2)
                    for j = 1:size(datcat(catk).segment_lengths{movk}{i},1)
                        cum_segment_lengths(j) = cum_segment_lengths(j) + datcat(catk).segment_lengths{movk}{i}(j);
                    end
                end
            end
            
            datcat(catk).inst_vel{movk} = datmovk.traj.inst_vel;
            datcat(catk).proc_vel{movk} = datmovk.traj.proc_vel;
            % datcat(catk).loc_alpha{movk} = datmovk.traj.loc_alpha;
            
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
        figure,totrl=gcf; %initialize figure
        [trl_n, trl_edges]=histcounts(datcat(catk).cum_run_length, 'BinWidth', rl_binwidth, 'Normalization', 'pdf');
        nhist_trl=trl_n;
        xhist_trl=trl_edges+(rl_binwidth/2);
        xhist_trl(end)=[]; 
        figure(totrl), hold on 
        bar(xhist_trl,nhist_trl)
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
        figure,totuvel=gcf; %initialize figure
        [uvel_n, uvel_edges]=histcounts(datcat(catk).cum_mean_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
        nhist_uvel=uvel_n;
        xhist_uvel=uvel_edges+(vel_binwidth/2);
        xhist_uvel(end)=[]; 
        figure(totuvel), hold on 
        bar(xhist_uvel,nhist_uvel)
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
        figure,boundtime=gcf; %initialize figure
        [bt_n, bt_edges]=histcounts(datcat(catk).cum_association_time, 'BinWidth', time_binwidth, 'Normalization', 'pdf');
        nhist_bt=bt_n;
        xhist_bt=bt_edges+(time_binwidth/2);
        xhist_bt(end)=[]; 
        figure(boundtime), hold on 
        bar(xhist_bt,nhist_bt)
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
        figure,localalpha=gcf; %initialize figure
        [loca_n, loca_edges]=histcounts(datcat(catk).cum_loc_alpha, 'BinWidth', loca_binwidth, 'Normalization', 'pdf');
        nhist_loca=loca_n;
        xhist_loca=loca_edges+(loca_binwidth/2);
        xhist_loca(end)=[]; 
        opt = statset('MaxIter',1e5,'Display','iter','TolFun',1e-20);
        sigmastart = cat(3, 0.01, 0.4);
        p0 = struct('mu', [0; 2], 'Sigma', sigmastart, 'ComponentProportion', [0.5; 0.5]);
        gmm_loc_alpha = fitgmdist(datcat(catk).cum_loc_alpha, 2, 'start', p0, 'options', opt); %
        y_loca = pdf(gmm_loc_alpha,xhist_loca');
        figure(localalpha), hold on 
        bar(xhist_loca,nhist_loca)
        plot(xhist_loca,y_loca,'-')
        xlabel('Local alpha-value'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Local alpha-value histogram'])
        hold off
        
        %instantaneous velocity
        figure,instvel=gcf; %initialize figure
        [instvel_n, instvel_edges]=histcounts(datcat(catk).cum_inst_vel, 'BinWidth', vel_binwidth, 'Normalization', 'pdf');
        nhist_instvel=instvel_n;
        xhist_instvel=instvel_edges+(vel_binwidth/2);
        xhist_instvel(end)=[]; 
        opt = statset('MaxIter',1e5,'Display','iter','TolFun',1e-20);
        sigmastart = cat(3, 0.1, 0.1);
        p0 = struct('mu', [0; 1500], 'Sigma', sigmastart, 'ComponentProportion', [0.2; 0.8]);
        gmm_inst_vel = fitgmdist(datcat(catk).cum_inst_vel', 2, 'options', opt); %'start', p0, 
        y_instvel = pdf(gmm_inst_vel,xhist_instvel');
        figure(instvel), hold on 
        bar(xhist_instvel,nhist_instvel)
        plot(xhist_instvel,y_instvel,'-')
        xlabel('Instantaneous velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Instantaneous velocity histogram'])
        hold off
        
        %instantaneous processive velocity
        figure,procvel=gcf; %initialize figure
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
        figure(procvel), hold on 
        bar(xhist_procvel,nhist_procvel)
        xlabel('Processive velocity (nm/s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Processive velocity histogram'])
        yprocV = normpdf(xhist_procvel, estimprocV(1), estimprocV(2));
        yprocVlo = normpdf(xhist_procvel, pciprocV(1,1), pciprocV(1,2));
        yprocVhi = normpdf(xhist_procvel, pciprocV(2,1), pciprocV(2,2));
        plot(xhist_procvel,1.0.*yprocV,'r-');
        plot(xhist_procvel,1.0.*yprocVlo,'r.');
        plot(xhist_procvel,1.0.*yprocVhi,'r.');
        hold off
        
        % landing position along MT - distance from plus-end
        figure,landpos = gcf;
        [landpos_n, landpos_edges]=histcounts(datcat(catk).cum_landing_dist_to_mt_end, 'BinWidth', 300, 'Normalization', 'pdf');
        nhist_landpos=landpos_n;
        xhist_landpos=landpos_edges+(500/2);
        xhist_landpos(end)=[]; 
        figure(landpos), hold on 
        bar(xhist_landpos,nhist_landpos)
        xlabel('Landing distance to MT plus-end (nm)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Landing distance from MT plus-end'])
        hold off
        
        % landing position along MT - normalized
        figure,normlandpos = gcf;
        [normlandpos_n, normlandpos_edges]=histcounts(datcat(catk).cum_norm_landing_pos, 'BinWidth', 0.1, 'Normalization', 'pdf');
        nhist_normlandpos=normlandpos_n;
        xhist_normlandpos=normlandpos_edges+(0.1/2);
        xhist_normlandpos(end)=[]; 
        figure(normlandpos), hold on 
        bar(xhist_normlandpos,nhist_normlandpos)
        xlabel('Landing position along MT (fraction of 1)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Normalized Landing Position'])
        hold off
        
        % landing time on given MT
        figure, landtime = gcf;
        figure(landtime), hold on
        mt_lengths = [];
        num_mov = size(datcat(catk).track_start_times,2);
        mtnum = 0;
        for movj = 1:num_mov
            num_mt_mov = size(datcat(catk).track_start_times{1,movj},1);
            for mtj = 1:num_mt_mov
                mtnum = mtnum+1;
                mt_length_temp = (datcat(catk).mt_lengths{movj}{mtj})/1000;
                if ~isempty(mt_length_temp)
                    mt_lengths(mtnum) = mt_length_temp;
                end
                num_tracks_mt = size(cell2mat(datcat(catk).track_start_times{1,movj}(mtj,1)),1);
                cmap = colormap(parula(num_mov));
                scatter(repmat(mtnum,[num_tracks_mt,1]), cell2mat(datcat(catk).track_start_times{1,movj}(mtj,1)).*mt_lengths(mtnum),25,cmap(movj,:),'filled')
            end 
        end
        xlabel('MT number'), ylabel('Track start time (s*\mum)'), title([motor{mk},' ', mt_type{mtk},' Landing times on each MT'])
%         xlim([0 mtnum])
        hold off
        
        figure, timebwland = gcf;
        figure(timebwland), hold on
        num_mov = size(datcat(catk).track_start_times,2);
        mtnum = 0;
        for movj = 1:num_mov
            num_mt_mov = size(datcat(catk).track_start_times{1,movj},1);
            for mtj = 1:num_mt_mov
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
                    xlabel('MT number'), ylabel('Time between track start events (s*\mum)'), title([motor{mk},' ', mt_type{mtk},' Time between start of tracks'])
                    xlim([0 mtnum])
                    %hold off
                end
            end 
        end
        figure,timebwlandhist = gcf;
        [deltland_n, deltland_edges]=histcounts(cum_time_bw_landing, 'BinWidth', 0.5, 'Normalization', 'pdf');
        nhist_deltland=deltland_n;
        xhist_deltland=deltland_edges+(0.5/2); %%%%%%%%%%%%%%%
        xhist_deltland(end)=[]; 
        figure(timebwlandhist), hold on 
        bar(xhist_deltland,nhist_deltland)
        xlabel('Time between landing events (s)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Time between landing events'])
        hold off
        
        figure,timebwlandbydisthist = gcf;
        [deltlanddist_n, deltlanddist_edges]=histcounts(cum_time_bw_landing_by_length, 'BinWidth', 1, 'Normalization', 'pdf');
        nhist_deltlanddist=deltlanddist_n;
        xhist_deltlanddist=deltlanddist_edges+(1/2);
        xhist_deltlanddist(end)=[]; 
        figure(timebwlandbydisthist), hold on 
        bar(xhist_deltlanddist,nhist_deltlanddist)
        xlabel('Time between landing events (s*\mum)'), ylabel('Probability density'), title([motor{mk},' ', mt_type{mtk},' Time between landing events'])
        hold off
        
        %track start segment
        if zcap == 1
            figure,trackstartsegment=gcf; %initialize figure
            [tkstartseg_n, tkstartseg_edges]=histcounts(datcat(catk).cum_track_start_segment, 'BinEdges', [0.5,1.5,2.5,3.5,4.5,5.5], 'Normalization', 'count');
            nhist_tkstartseg=tkstartseg_n;
            xhist_tkstartseg= [1,2,3,4,5];
            %xhist_tkstartseg(end)=[]; 
            for j = 1:size(xhist_tkstartseg,2)
                nhist_tkstartseg(j) = nhist_tkstartseg(j)/(cum_segment_lengths(j)/1000);
            end
            figure(trackstartsegment), hold on 
            bar(xhist_tkstartseg,nhist_tkstartseg)
            set(gca,'XTick',[1,2,3,4,5], 'xticklabel',{'Plus GMP-CPP cap','Plus GDP lattice','GMP-CPP seed','Minus GDP lattice','Minus GMP-CPP cap'});
            xlabel('MT Segment'), ylabel('Number of tracks starting in segment per unit distance (/\mum)'), title([motor{mk},' ', mt_type{mtk},' Track Start Segment'])

            figure, trackstartsegmenttype=gcf;
            nhist_tk_starttype = zeros(3,1);
            xhist_tkstarttype = [1,2,3];
            figure(trackstartsegmenttype), hold on 
            nhist_tkstarttype(1) = nhist_tkstartseg(1)+nhist_tkstartseg(5);
            nhist_tkstarttype(2) = nhist_tkstartseg(2)+nhist_tkstartseg(4);
            nhist_tkstarttype(3) = nhist_tkstartseg(3);
            bar(nhist_tkstarttype)
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