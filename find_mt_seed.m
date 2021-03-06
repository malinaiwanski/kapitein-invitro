%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski 2020-03-20
%% This function determines the orientation of the MT and which indices of the MT are part of the MT seed
% input: mt_coords (interpolated x,y coordinates of points along MT), seg_boundaries (x,y coordinates ON MT of the segment boundaries), and response (output by Curve Tracing plugin for Fiji available here: https://github.com/jalmar/CurveTracing)
% output: flip_mt (1 if MT coordinates as listed are from + to - end; 0 if MT coordinates as listed are from - to + end) and mt_seed (indices of MT defined as being part of the seed -based on a more brightly labelled seed))

function [flip_mt, mt_seed] = find_mt_seed(mt_coords,seg_boundaries,response)
    tot_mt_length = arclength(mt_coords(:,1),mt_coords(:,2));
    
    %plot([1:1:length(response)],response)
    
    findchangepts(response,'MinThreshold',0.3)
    
    
    
    %order segment boundaries based on distance from plus-end of MT
    if ~isempty(seg_boundaries)
        seg_boundaries_old = boundaries_on_mt;
        dist_to_bound = double.empty(length(boundaries_on_mt),2);
        for j=1:size(seg_boundaries,1)
            [~,interp_mt_bound_ind] = ismember(seg_boundaries(j,:),mt_coords, 'rows');
            interp_mt_bound_ind = nonzeros(interp_mt_bound_ind);
            mt_inds_to_bound = interp_mt_bound_ind:1:size(mt_coords,1);
            dist_to_bound(j,1) = arclength(mt_coords(mt_inds_to_bound,1),mt_coords(mt_inds_to_bound,2)); %distance to plus-end
            dist_to_bound(j,2) = j;
        end
        dist_to_bound = sortrows(dist_to_bound,1);
        for j=1:size(seg_boundaries,1)
            seg_boundaries(j,:) = seg_boundaries_old(dist_to_bound(j,2),:); %boundaries ordered from closest to furthest from plus-end of MT 
        end
        dist_along_mt = [0;dist_to_bound(:,1);tot_mt_length];

        % find length of each segment
        segment_lengths = abs(diff(dist_along_mt)); %either 3 or 5 entries in order +cap, +gdp, seed(, -gdp, -cap)
    end
    
    %%%% temp
    flip_mt = 0;
    mt_seed = [1,2,3,4];
    
%     if ~isempty(seg_boundaries)
%         seg_boundaries_old = seg_boundaries;
%         for j=1:size(seg_boundaries,1)
%             [~,interp_mt_bound_ind] = ismember(seg_boundaries(j,:),mt_coords, 'rows');
%             interp_mt_bound_ind = nonzeros(interp_mt_bound_ind);
%             mt_inds_to_bound = 1:1:interp_mt_bound_ind; %:1:size(mt_coords,1);
%             dist_to_bound(j,1) = arclength(mt_coords(mt_inds_to_bound,1),mt_coords(mt_inds_to_bound,2)); %distance from start of MT
%             dist_to_bound(j,2) = j;
%         end
%         dist_to_bound = sortrows(dist_to_bound,1);
%         for j=1:size(seg_boundaries,1)
%             segment_intensity = mean(response(
%         end
%     end
%     dist_along_mt = [0;dist_to_bound(:,1);tot_mt_length];
% 
%     % find length of each segment
%     segment_lengths{mttk} = abs(diff(dist_along_mt)); %either 3 or 5 entries in order +cap, +gdp, seed(, -gdp, -cap)
% 
%     if size(dist_along_mt,1) == 4 
%         segment_annotate = 3;
%     elseif size(dist_along_mt,1) == 6
%         segment_annotate = 5;
%     else
%         segment_annotate = 0;
%         disp(['MT number ', num2str(mttk), ' does not have 3 or 5 segments. Retrace MT and/or segment boundaries']);
%     end

    
end