%% UU - Kapitein Lab
% Analyze in vitro single molecule motility assays
% MK Iwanski and C Chen 2019-11-05

%% Original code: McGill - Hendricks Lab
% MK Iwanski 2019-05-29
% with AG Hendricks 

%% This function calculates the MSD of a trajectory in 2D
% inputs:
    % r is a cell array of position data, each entry in r is a separate trajectory (i.e. r{1}=[x1(1),y1(1); x1(2),y1(2); x1(3),y1(3); ...]; r{2}=[x2(1),y2(1); x2(2),y2(2); x2(3),y2(3); ...];)
    % delays is an array of delays (in units of the index of r) - delay must be less than the number of points in the trajectory
% outputs:
    % y is MSD for each delta_t
    % y_std is standard deviation for each delta_t
    % N is number of points for each delta_t
    % msd_delay1 is MSD for all delta_t = 1

function [y,y_std,N,msd_delay1]=MSD_2D(r,delays)

%%
msd_delay1 = [];
cr=1; %[pixel-size correction]
Ndat=numel(r); %no. of trajectories
K=delays;%1:1:10;

for j=1:numel(K) %for each delay
    k=K(j);
    dispsqk=[];
    for kd=1:Ndat %for each trajectory
        clear rk
        rk = [];
        rk=r{kd}; 

        [mr,nr]=size(rk); %mr will be the number of points in the trajectory; nr should always be 2 (x and y)

        if mr<nr %if input had rows and columns switched, switch them for calculation
            rk=rk'; 
            [mr,nr]=size(rk);
        end
        if k<(mr-1)
            jmax=min(mr,max(K)*4);
            j0=1:(jmax-k);
            j1=(k+1):jmax;
            dkd = (rk(j1,1)-rk(j0,1)).^2+(rk(j1,2)-rk(j0,2)).^2; %calculate MSD
            dispsqk=[dispsqk;dkd]; 
        end
        if k == 1
            msd_delay1 = [msd_delay1;dispsqk]; %msd for 1 frame delay
        end
    end
    jk=find(isnan(dispsqk)==0);
    N(j)=numel(dispsqk(jk)); %no. of sample points for each delta_t
    y(j)=sum(dispsqk(jk))/N(j); %mean MSD for each delta_t
    y_std(j)=std(dispsqk(jk)); %standard deviation for each delta_t
end


end