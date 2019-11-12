%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski and C Chen 2019-11-05

%% Original code: McGill - Hendricks Lab
% MK Iwanski 2019-05-29
% with AG Hendricks (transient MSD over sliding window 2019-04-10)

%% This function finds local alpha-values by calculating MSD over a sliding window
% It assigns the calculated alpha-value to the middle point in the window (except at end points the first and last calculated alpha-values are repeated)
% Input: 
    % x (x positions of trajectory [nm])
    % y (y positions of trajectory [nm])
    % frames (frame of each localization in trajectory)
    % l_window (number of frames over which MSD is locally calculated)
    % exp_time (time between frames [s])
    % msd_thresh (alpha-value above which is processive, below which is paused)
    % msd_step (minimum threshold for findchangepts function; minimum improvement in residual error; changepoints currently based on mean values, but can use rms, std, mean and slope)
    % l_min (smallest no. of frames a section can be; minimum distance between two changepoints)
% Output: 
    % [frame number, local alpha-value, processive (1) or paused (0)]

function tmsd_changepts = tMSD_2D(x,y,frames,l_window,exp_time,msd_thresh, msd_step,l_min)
r = [x, y];          
t = frames.*exp_time;

% vels = sqrt(diff(x).^2+diff(y).^2)./diff(t); %calculate frame-to-frame velocity using interpolated data

delays = floor(l_window/4):ceil(l_window/2); %calculate which delays to use based on size of averaging window
logdelays = log10(delays.*exp_time);

Nk = numel(r(:,1)) - l_window; %number of windows
mmsd = zeros(numel(t),2); %to store local alpha-values

for k=1:Nk
    jk = k:(k+l_window);
    MSDk = MSD_2D({r(jk,:)},delays); %calculate MSD of trajectory within window
    pk = polyfit(logdelays,log10(MSDk),1); %calculate alpha-value using linear fit to log-log plot of MSD over time
    mmsd(k+ceil(l_window/2),:) = pk; %slope is alpha-value (first column)
end

for i = 1:ceil(l_window/2)
    mmsd(i,:) = mmsd(ceil(l_window/2)+1,:); %assign alpha-value to middle point in window
end
for i = size(mmsd,1)-l_window+ceil(l_window/2)+1:size(mmsd,1)
    mmsd(i,:) = mmsd(size(mmsd,1)-l_window+ceil(l_window/2),:); %for end points, repeat initial or final calculated alpha-value
end

%% plot displacement of trajectory and corresponding alpha-values before parsing
% figure
% disp = sqrt(sum((r-r(1,:)).^2,2)); 
% subplot(2,1,1), plot(t,disp), ylabel('Displacement')
% subplot(2,1,2), plot(t,mmsd(:,1)), xlabel('Time (s)'), ylabel('Alpha')

%% findchangepts 
% %plots of means and residuals
% figure
% findchangepts(mmsd(:,1),'MinThreshold',msd_step,'MinDistance',l_min)

proc = zeros(size(mmsd,1),1);
% changepts = findchangepts(abs(mmsd(:,1)),'MinThreshold',msd_step,'MinDistance',l_min); %version using absolute value
changepts = findchangepts((mmsd(:,1)),'MinThreshold',msd_step,'MinDistance',l_min); %changepoints identified are only used as putative places where transition paused <--> processive
num_changes = length(changepts);
change_ind = unique([1;changepts;size(mmsd,1)]); %add first and last frame
num_changes = length(change_ind);

%% parsing
for i =1:(num_changes)
    if i == num_changes %for last section
%         meanmsd = mean(abs(mmsd(change_ind(i):1:change_ind(end),1))); %version using absolute value
%         msdquantile = quantile(abs(mmsd(change_ind(i):1:change_ind(end))),[0.01 0.95]); %version using absolute values
        
        meanmsd = mean((mmsd(change_ind(i):1:change_ind(end),1))); %mean alpha value between last changepoint and last frame
        
        msdquantile = quantile((mmsd(change_ind(i):1:change_ind(end))),[0.01 0.95]); %not currently used
        msdspan = msdquantile(2)-msdquantile(1); %not currently used

        if meanmsd >= msd_thresh % if mean alpha between sequential changepoints is above threshold --> processive
            proc(change_ind(i):1:change_ind(end),1) = 1;
        end
        
    else %for all other sections
%         meanmsd = mean(abs(mmsd(change_ind(i):1:change_ind(i+1),1))); %version using absolute value
%         msdquantile = quantile(abs(mmsd(change_ind(i):1:change_ind(i+1))),[0.01 0.95]); %version using absolute values
        meanmsd = mean((mmsd(change_ind(i):1:change_ind(i+1),1))); %mean alpha value between sequential changepoints
        
        msdquantile = quantile((mmsd(change_ind(i):1:change_ind(i+1))),[0.01 0.95]); %not currently used
        msdspan = msdquantile(2)-msdquantile(1); %not currently used
    %     if ((velspan >= 0.09)) %alternate analysis strategy
    %         proc(change_ind(i):1:change_ind(i+1),1) = 1;
    %     end
        if meanmsd >= msd_thresh % if mean alpha between sequential changepoints is above threshold --> processive
            proc(change_ind(i):1:change_ind(i+1),1) = 1;
        end
    end
end
proc_frames = find(proc);
pause_frames = find(~proc);
frames = [proc_frames;pause_frames];
frames = sort(frames);

%%
tmsd_changepts = [frames,mmsd(:,1),proc];

% %% plot parsed trajectory: displacement over time and local alpha-values over time; paused in red, processive in blue
% figure 
% disp = sqrt(sum((r-r(1,:)).^2,2)); 
% subplot(2,1,1), hold on
% plot(proc_frames,disp(proc_frames),'b.') 
% plot(pause_frames,disp(pause_frames),'r.') 
% ylabel('Displacement')
% hold off
% subplot(2,1,2), hold on
% plot(proc_frames,mmsd(proc_frames,1),'b.')
% plot(pause_frames,mmsd(pause_frames,1),'r.')
% xlabel('Time (s)'), ylabel('Alpha')
% hold off

end
