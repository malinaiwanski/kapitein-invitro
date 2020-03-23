%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski 2020-03-20
%% This function identifies the mean velocity of different types of motility and/or on different MT lattices
% input: motorqpos (position of motor ON MT), position (x,y localizations
% of motor), frames (frames of motor localizations in movie), interp_mt
% (interpolated x,y MT coordinates), seg_ind (vector identifying in which
% frames motor is in a condition (e.g. cap and processive) with a 1,
% otherwise a 0)
% output: vector containing mean velocity of each segment (of type
% determined by the input vector seg_ind (e.g. if seg_ind specifies
% processive frames on the GDP lattice, output mean velocities will be for
% processive segments of the trajectory on the GDP lattice)

function [mu_segvels] = mean_segment_vel(motorqpos, position, frames, exp_time, interp_mt, seg_ind)
    
    %intialize variables
    motorq_mtind = [];
    mu_segvels = [];                       
    seg_vel = [];

    changes = diff(seg_ind); %-1 means changes out of state; 1 means changes into state
    changeframes =find(changes); %gives the last frame of each segment - these are not absolute frames, but the first frame of a track is always 1
    segment_ends = [changeframes;size(seg_ind,1)]; %the last frame of each segment including the last frame - these are not absolute frames, but the first frame of a track is always 1
    segment_starts = [1; changeframes+1];

    ind_to_check = sort([segment_ends;segment_starts]);
    for ik = 1:size(ind_to_check,1)
        [~,motor_mtind] = ismember(interp_mt, motorqpos(ind_to_check(ik),:), 'rows');
        motor_mtind = find(motor_mtind);
        motorq_mtind = [motorq_mtind;motor_mtind];
    end

    pn = 0;
    rn = 0;
    if (isempty(changeframes) && ~isempty(seg_ind)) %motor does not switch between paused and processive
        if seg_ind(1) == 0
            mu_segvels = [];
        elseif seg_ind(1) == 1 
            seg_mt_ind = motorq_mtind(1):1:motorq_mtind(end);
            if length(seg_mt_ind) > 1
                seg_disp = arclength(interp_mt(seg_mt_ind,1),interp_mt(seg_mt_ind,2));
                usegvel = seg_disp/(size(frames,1)*exp_time);
                mu_segvels = [mu_segvels; usegvel];
                seg_vel = diff(position)./diff(frames.*exp_time);
                rn = rn+1;
            end
        end
    elseif(~isempty(changeframes)) %motor does switch
        spans = [changeframes(1); diff(changeframes); length(changes)+1-changeframes(end)];
        for sn = 1:(numel(spans))
            if sn == numel(spans)
                if seg_ind(end) == 0
                    segment_ind = segment_starts(sn):1:size(seg_ind,1);
                    pn = pn  + 1;
                elseif seg_ind(end) == 1 
                    seg_mt_ind = motorq_mtind(2*sn-1):1:motorq_mtind(2*sn);
                    if length(seg_mt_ind) > 1
                        segment_ind = segment_starts(sn):1:size(seg_ind,1);
                        seg_vel = [seg_vel, diff(position(segment_ind))./diff(frames(segment_ind)'.*exp_time)]; %changed from run_vel
                        seg_disp = arclength(interp_mt(seg_mt_ind,1),interp_mt(seg_mt_ind,2));
                        usegvel = seg_disp/(spans(sn)*exp_time);
                        mu_segvels = [mu_segvels; usegvel];
                        rn = rn + 1;
                    end
                end
            else
                if (changes(changeframes(sn)) == 1) %was paused
                    segment_ind = segment_starts(sn):1:segment_starts(sn+1)-1;
                    pn = pn + 1;
                elseif (changes(changeframes(sn)) == -1) %was processive
                    %store run time lengths
                    seg_mt_ind = motorq_mtind(2*sn-1):1:motorq_mtind(2*sn);
                    if length(seg_mt_ind) > 1
                        segment_ind = segment_starts(sn):1:segment_starts(sn+1)-1;
                        seg_vel = [seg_vel, diff(position(segment_ind))./diff(frames(segment_ind)'.*exp_time)];  %changed from run_vel                                        
                        seg_disp = arclength(interp_mt(seg_mt_ind,1),interp_mt(seg_mt_ind,2));
                        usegvel = seg_disp/(spans(sn)*exp_time);
                        mu_segvels = [mu_segvels; usegvel];
                        rn = rn + 1;
                    end
                end
            end
        end
    end
end