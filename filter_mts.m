%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MKI and CC 2019-11-05
%% This function takes in a cell of all MTs in a movie, formats, and filters them
% input: Cell with columns %1=id %2=roi_name %3=x %4=y, filtering options, pixel size
% output: Filtered cell where each mt{i} is an array with columns %1=x %2=y corresponding to one mt that meets the criteria (e.g. does not cross another mt, etc.)

function [temp_mts] = filter_mts(mt_data, filt_cross, cross_dist, pixel_size, zplot)
% initialize figure if plotting
if zplot ~= 0
    figure, mt_plot=gcf;
    xlabel('x (nm)'), ylabel('y (nm)'), title('Microtubules')
end

%format MT xy positions
mts_id = str2double(mt_data{1,1}(:))+1;
mts_x = str2double(mt_data{1,3}(:));
mts_y = str2double(mt_data{1,4}(:));

num_mts = length(unique(mts_id));

temp_mts = cell(num_mts,1);
for i = 1:num_mts
    skip_mts = zeros(1,num_mts); %vector to store which MTs do not meet filtering criteria (1) or do meet it (0)
    mt_id = find(mts_id(:,1)==i,1);
    mt_l = length(find(mts_id(:,1)==i));
    temp_mts{i} = [mts_x([mt_id:1:mt_id+mt_l-1],1).*pixel_size+2*pixel_size, ys([mt_id:1:mt_id+mt_l-1],1).*pixel_size+2*pixel_size]; %MT position x,y in um
    %add one pixel dimension because ImageJ starts at 0, MATLAB at 1; add another because positions are plotted at the bottom,left of a pixel (??)

    %FILTER: skip MTs with single unique x-value (vertical MTs)
    [~,uni_ind] = unique(temp_mts{i}(:,1),'stable'); %find repeated x-values in MT
    uni_ind = sort(uni_ind);
    if length(uni_ind) == 1 
        %uni_ind
        skip_mts(i) = 1;
        continue %do not need to filter
    end
    
    if filt_cross ~= 0 %filter out MTs that are near other MTs
        %interpolate MT coordinates
        x_interp = temp_mts{i}(1,1):0.01:temp_mts{i}(end,1);
        interp_mts{i} = [];
        interp_mts{i}(:,1)=interp1(temp_mts{i}(uni_ind,1),temp_mts{i}(uni_ind,1),x_interp,'linear'); %interpolated x coordinates of MT
        interp_mts{i}(:,2)=interp1(temp_mts{i}(uni_ind,1),temp_mts{i}(uni_ind,2),x_interp,'linear'); %interpolated y coordinates of MT
        interp_mtcoords = [interp_mtcoords;interp_mts{i}];

        %find points very close (defined by <= cross_dist) to other MTs (i.e. regions of MT overlap)
        mt_kq =[xk_correct,yk_correct];
        other_mtsq = setdiff(interp_mtcoords,interp_mts{i,1},'rows','stable');
        [points_near_diff_mt,~] = rangesearch(mt_kq,other_mtsq,cross_dist,'Distance','euclidean');
        points_near_diff_mt = points_near_diff_mt(~cellfun('isempty',points_near_diff_mt));
        for ip = 1:size(points_near_diff_mt,1)
            points_near_diff_mt{ip}=points_near_diff_mt{ip}';
        end
        points_near_diff_mt = cell2mat(points_near_diff_mt);
        points_near_diff_mt = unique(points_near_diff_mt,'stable');
        if ~isempty(points_near_diff_mt)
            mt_near_diff_mt = [mt_kq(points_near_diff_mt,1),mt_kq(points_near_diff_mt,2)];
            skip_mts(i) = 1;
        end
    end
   
    if zplot ~=0
        figure(mt_plot), hold on, plot(interp_mts{i}(:,1),interp_mts{i}(:,2),'-','Color',[0 0 0])
        figure(mt_plot), hold on, plot(mt_near_diff_mt(:,1),mt_near_diff_mt(:,2),'*','Color',[0 0 0]) %parts of MT considered crossing
    end
end
mts = temp_mts;
mts{skip_mts} = {};
mts = mts(~cellfun('isempty',mts));

end