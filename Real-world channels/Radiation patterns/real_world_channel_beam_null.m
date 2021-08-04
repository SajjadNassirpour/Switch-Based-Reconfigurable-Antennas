% -- Generate the real-world complex channels based on ray-tracing model -- %

% K :                   Number of users
% M :                   Number of elements in array antenna
% iters :               Number of trials
% est_SNR :             SNR(dB) to estimate the channels
% user_distr :          Distribution of users around the BTSs

function [kbyk_chans,kbyk_chans_nonoise,kbyk_analogbeamformvecs,kbyk_analogbeamformvecs_nulls,kbyk_analogbeamformvecs_nulls_multipath,bts_cords,all_selected_cords,grid_dim,bts_cord_ants,refl_cords]=real_world_channel_beam_null(K,M,iters,est_SNR,user_distr)

% Number of channels
num_chans = iters;

% Generate sample points
geom_spacing = 100;
grid_dim = 15;
sample_points_vec = linspace(0.5,grid_dim-0.5,geom_spacing);
[X,Y] = meshgrid(sample_points_vec,sample_points_vec);
sample_points = [X(:) Y(:)];

%Num users
num_users = K;
num_choices = geom_spacing^2;

% BTS cords
num_bts = num_users; % k*k channel
x_offset = grid_dim/4;
y_cord1 = grid_dim/4;
y_cord2 = (grid_dim*3)/4;
x_sep = grid_dim/(num_users/2);
bts_cords = ones([num_users,2]);
bts_cords(1:1:num_users/2,2) = y_cord1;
bts_cords(1:1:num_users/2,1) = x_offset+(0:1:(num_users/2-1))*x_sep;
bts_cords(num_users/2+1:1:num_users,2) = y_cord2;
bts_cords(num_users/2+1:1:num_users,1) = x_offset+(0:1:(num_users/2-1))*x_sep;


% num_ants
num_ants = M;

bts_cord_ants = zeros(num_chans,num_bts,num_ants,2);

kbyk_chans = ones(num_chans,num_bts,num_users,num_ants);
kbyk_chans_nonoise = ones(num_chans,num_bts,num_users,num_ants);

kbyk_analogbeamformvecs = ones(num_chans,num_bts,num_ants);
kbyk_analogbeamformvecs_nulls = ones(num_chans,num_bts,num_ants);

kbyk_analogbeamformvecs_nulls_multipath = ones(num_chans,num_bts,num_ants);

chan_est_snr_db = est_SNR;

dist_thresh = 0.5;

chan_index=1;

all_selected_cords = zeros(num_chans,num_users,2);
refl_cords = zeros(num_chans,num_users,num_users,4,2);

user_distr ='spatial'; % can be spatial or uniform

while(chan_index<=num_chans)
    
    if(user_distr == 'uniform')
        selected_user_cords = sample_points(randperm(num_choices,num_users),:);
    elseif(user_distr == 'spatial')
        selected_user_cords = zeros(num_users,2);
        for k=1:1:num_users
            y_grid_ind = floor((k-1)/2);
            x_grid_ind = mod(k-1,2);
            
            x_grid_min = max((x_grid_ind)*grid_dim/4+0.5-grid_dim/4,0);
            x_grid_max = min((x_grid_ind+1)*grid_dim/4-0.5+grid_dim/4,grid_dim);
%             
%             x_grid_min = 0;
%             x_grid_max = grid_dim;
            
            y_grid_min = (y_grid_ind)*grid_dim/2+0.5;
            y_grid_max = (y_grid_ind+1)*grid_dim/2-0.5;
            y_grid_min = (y_grid_min+y_grid_max)/2;
            
            selected_user_cords(k,:) = [unifrnd(x_grid_min,x_grid_max), unifrnd(y_grid_min,y_grid_max)];
            
        end
    end
    % calculate distance to each ap, remove cases where coordinates were
    % very close
    unselect = 0;
    for j=1:1:num_bts
        curr_bts_cords = bts_cords(j,:);
        curr_dists = vecnorm(selected_user_cords-curr_bts_cords,2,2);
        min_curr_dist = min(curr_dists);
        if(min_curr_dist<dist_thresh)
            unselect =1;
            break
        end
    end
    if(unselect==1)
        continue
    end
        
    manual_cords= [[7,3.5];[14.5,10.8];[10,4.5];[13.5,11.5];]; % almost
%     manual_cords = [[0.3,4.5];[6.2,3.3];[0.5,14];[7,11.8]];
    selected_user_cords = manual_cords;
    all_selected_cords(chan_index,:,:) = selected_user_cords;
%     all_selected_cords(chan_index,:,:) = selected_user_cords;

    for j=1:1:num_bts
        curr_bts_cords = bts_cords(j,:);
        [kbyk_chans(chan_index,j,:,:), kbyk_chans_nonoise(chan_index,j,:,:), kbyk_analogbeamformvecs(chan_index,j,:,:), kbyk_analogbeamformvecs_nulls(chan_index,j,:,:),kbyk_analogbeamformvecs_nulls_multipath(chan_index,j,:,:),bts_cord_ants(chan_index,j,:,:),refl_cords(chan_index,j,:,:,:)] =  chan_est_geometrical_kbyk_null(selected_user_cords, curr_bts_cords, chan_est_snr_db, grid_dim,  num_ants, j);
%         kbyk_chans(chan_index,j,:,:)
    end
    
    chan_index= chan_index+1;
end

end

