function [mu_segvels] = mean_segment_vel(seg_ind)

    motorq_mtind = [];
    
    mean_run_vel = [];                       
    proc_vel = []; %inst_vel(proc_frames==1); %diff(position(proc_frames==1))./diff(frame_tk(proc_frames==1).*exp_time);
    pause_vel = []; %inst_vel(proc_frames==0); %diff(position(proc_frames==0))./diff(frame_tk(proc_frames==0).*exp_time);

    changes = diff(seg_ind); %-1 means changes processive to paused; 1 means changes paused to processive
    changeframes =find(changes); %gives the last frame of each segment - these are not absolute frames, but the first frame of a track is always 1
    segment_ends = [changeframes;size(seg_ind,1)]; %the last frame of each segment including the last frame - these are not absolute frames, but the first frame of a track is always 1
    segment_starts = [1; c hangeframes+1];

    ind_to_check = sort([segment_ends;segment_starts]);
    motorqpos = motor_on_mt{ftk}; %%%%%
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