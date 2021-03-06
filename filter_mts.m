%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski 2019-11-05
%% This function takes in a cell of all MTs in a movie, formats, and filters them
% input: Cell with columns %1=id %2=roi_name %3=x %4=y, filtering options, pixel size
% output: Filtered cell where each mt{i} is an array with columns %1=x %2=y corresponding to one mt that meets the criteria (e.g. does not cross another mt, etc.)

function [mts, interp_mts, response_mts, skip_mts] = filter_mts(mt_data, analyze_mt_num, filt_cross, cross_dist, filt_short, min_length, pixel_size, num_pix_x, num_pix_y, zplot, zcap, zcurve)

% initialize figure if plotting
if zplot ~= 0
    figure, mt_plot=gcf;
    xlabel('x (nm)'), ylabel('y (nm)'), title('Microtubules')
end

all_mt_near_diff_mt = [];

if analyze_mt_num ~= -1 %only analyze specified MT
    disp(strcat('Analyzing only MT ',num2str(analyze_mt_num)))
    filt_cross = 0; %if analyzing only one MT, do not filter out e.g. crossing MTs
    
    %format MT xy positions
    mts_id = str2double(mt_data{1,1}(:));
    if zcurve == 0
        mts_id = mts_id +1;
        mt_id = find(mts_id(:,1)==analyze_mt_num+1,1); %analyze_mt_num specified as starting from 0, but here IDs start at 1
    else
        mt_id = find(mts_id(:,1)==analyze_mt_num,1); %analyze_mt_num specified as starting from 0, but here IDs start at 1
        if zcap == 1
            mt_response = str2double(mt_data{1,5}(:));
        end
    end
    
    mts_x = str2double(mt_data{1,3}(:));
    mts_y = str2double(mt_data{1,4}(:));
    if zcap ==1
        mts_x = str2double(mt_data{1,4}(:));
        mts_y = str2double(mt_data{1,3}(:));
    end
    
    num_mts = 1;
else %analyze all MTs
    %format MT xy positions
    mts_id = str2double(mt_data{1,1}(:));
    if zcurve == 0
        mts_id = mts_id +1;
        mt_id = find(mts_id(:,1)==analyze_mt_num+1,1); %analyze_mt_num specified as starting from 0, but here IDs start at 1
    else
        mt_id = find(mts_id(:,1)==analyze_mt_num,1); %analyze_mt_num specified as starting from 0, but here IDs start at 1
        if zcap ==1
            mt_response = str2double(mt_data{1,5}(:));
        end
    end
    mts_x = str2double(mt_data{1,3}(:));
    mts_y = str2double(mt_data{1,4}(:));
    if zcap ==1
        mts_x = str2double(mt_data{1,4}(:));
        mts_y = str2double(mt_data{1,3}(:));
    end

    num_mts = length(unique(mts_id));
end
cmap=colormap(colorcube(num_mts));

% initialize variables
temp_mts = cell(num_mts,1);
response_mts = cell(num_mts,1);
interp_mts = cell(num_mts,1);
interp_mtcoords = [];
mt_near_diff_mt = cell(num_mts,1);

for i = 1:num_mts
    skip_mts = zeros(1,num_mts); %vector to store which MTs do not meet filtering criteria (1) or do meet it (0)
    cross_mts = zeros(1,num_mts);
    mt_id = find(mts_id(:,1)==i,1);
    mt_l = length(find(mts_id(:,1)==i)); %number of points making up MT, NOT physical length
    if zcurve == 0
        temp_mts{i} = [mts_x([mt_id:1:mt_id+mt_l-1],1).*pixel_size+1.5*pixel_size, mts_y([mt_id:1:mt_id+mt_l-1],1).*pixel_size+1.5*pixel_size]; %MT position x,y in nm
        %add one pixel dimension because ImageJ starts at 0, MATLAB at 1; add 0.5 because positions are plotted at the bottom,left of a pixel (??)
    else
        temp_mts{i} = [mts_x([mt_id:1:mt_id+mt_l-1],1).*pixel_size, mts_y([mt_id:1:mt_id+mt_l-1],1).*pixel_size]; %MT position x,y in nm
    end
    if temp_mts{i}(1,1) > temp_mts{i}(end,1)
        temp_mts{i} = flipud(temp_mts{i});
    end
    %temp_mts{i} = sortrows(temp_mts{i}); %put x-values in ascending order for later interpolation
    
    if zcap == 1 && zcurve == 1
        response_mts{i} = mt_response([mt_id:1:mt_id+mt_l-1],1);
    end
    
    %FILTER: skip MTs with single unique x-value (vertical MTs)
    [~,uni_ind] = unique(temp_mts{i}(:,1),'stable'); %find repeated x-values in MT
    %uni_ind = sort(uni_ind);
    if length(uni_ind) == 1 
        %uni_ind
        skip_mts(i) = 1;
    end
    
    %interpolate MT coordinates
    x_interp = temp_mts{i}(1,1):10:temp_mts{i}(end,1);
    interp_mts{i} = [];
    if skip_mts(i) == 0
        interp_mts{i}(:,1)=interp1(temp_mts{i}(uni_ind,1),temp_mts{i}(uni_ind,1),x_interp,'linear','extrap'); %interpolated x coordinates of MT
        interp_mts{i}(:,2)=interp1(temp_mts{i}(uni_ind,1),temp_mts{i}(uni_ind,2),x_interp,'linear','extrap'); %interpolated y coordinates of MT
    else
        interp_mts{i} = temp_mts{i};
    end
    interp_mtcoords = [interp_mtcoords;interp_mts{i}]; %accumulates all MT coordinates
end

%FILTER: skip MTs that are too short
if filt_short ~= 0
    for i = 1:num_mts
        %mt_length = sum(sqrt(temp_mts{i}(:,1).^2 + temp_mts{i}(:,2).^2));
        mt_length = arclength(temp_mts{i}(:,1),temp_mts{i}(:,2));
        if mt_length < min_length
            cross_mts(i) = 1;
        end
        if zcap == 1 && mt_length < 3*min_length
            cross_mts(i) = 1;
        end
    end
end

%FILTER: skip MTs that are near other MTs
if filt_cross ~= 0 
    for i = 1:num_mts
    %find points very close (defined by <= cross_dist) to other MTs (i.e. regions of MT overlap)
        mt_kq =[interp_mts{i}(:,1),interp_mts{i}(:,2)];
        other_mtsq = setdiff(interp_mtcoords,interp_mts{i,1},'rows','stable'); %compare to all MTs except for the current one
        [points_near_diff_mt,~] = rangesearch(mt_kq,other_mtsq,cross_dist,'Distance','euclidean');
        points_near_diff_mt = points_near_diff_mt(~cellfun('isempty',points_near_diff_mt));
        for ip = 1:size(points_near_diff_mt,1)
            points_near_diff_mt{ip}=points_near_diff_mt{ip}';
        end
        points_near_diff_mt = cell2mat(points_near_diff_mt);
        points_near_diff_mt = unique(points_near_diff_mt,'stable');
        if ~isempty(points_near_diff_mt)
            mt_near_diff_mt{i} = [mt_kq(points_near_diff_mt,1),mt_kq(points_near_diff_mt,2)];
            all_mt_near_diff_mt = [all_mt_near_diff_mt; mt_near_diff_mt{i}];
            %skip_mts(i) = 1;
            cross_mts(i) = 1;
        end
    end
end
  
if zplot ~=0
    for i = 1:num_mts
        figure(mt_plot), hold on, plot(interp_mts{i}(:,1),interp_mts{i}(:,2),'-','Color',cmap(i,:))
        if filt_cross ~= 0 && ~isempty(mt_near_diff_mt{i}) == 1 %~cellfun('isempty',mt_near_diff_mt{i})
            figure(mt_plot), hold on, plot(mt_near_diff_mt{i}(:,1),mt_near_diff_mt{i}(:,2),'*','Color',cmap(i,:)) %parts of MT considered crossing
        end
    end
end
    
% mts_to_skip = find(skip_mts == 1);
% mts = temp_mts;
% for j = 1:length(mts_to_skip)
%     mts{mts_to_skip(j)} = {};
%     interp_mts{mts_to_skip(j)} = {};
% end
% mts = mts(~cellfun('isempty',mts));
% interp_mts = interp_mts(~cellfun('isempty',interp_mts));
% cross_mts(mts_to_skip) = [];

mts = temp_mts;
skip_mts = skip_mts + cross_mts;
skip_mts = skip_mts~=0;
% cross_mts = skip_mts;

end