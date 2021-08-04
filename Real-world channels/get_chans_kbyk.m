rng(0)

% Number of channels
num_chans = 10;

% Generate sample points
geom_spacing = 100;
grid_dim = 15;
sample_points_vec = linspace(0.5,grid_dim-0.5,geom_spacing);
[X,Y] = meshgrid(sample_points_vec,sample_points_vec);
sample_points = [X(:) Y(:)];

%Num users
num_users = 8;
num_choices = geom_spacing^2;

% BTS cords
num_bts = num_users; % k*k channel
x_offset = grid_dim/8;
y_cord1 = grid_dim/4;
y_cord2 = (grid_dim*3)/4;
x_sep = grid_dim/(num_users/2);
bts_cords = ones([num_users,2]);
bts_cords(1:1:num_users/2,2) = y_cord1;
bts_cords(1:1:num_users/2,1) = x_offset+(0:1:(num_users/2-1))*x_sep;
bts_cords(num_users/2+1:1:num_users,2) = y_cord2;
bts_cords(num_users/2+1:1:num_users,1) = x_offset+(0:1:(num_users/2-1))*x_sep;


figure(1)
h=[];

h(1) = scatter(bts_cords(:,1),bts_cords(:,2),100,'s','filled');
hold on

ylim([-2,grid_dim+2]);
xlim([-2,grid_dim+2]);

% num_ants
num_ants = 16;
kbyk_chans = ones(num_chans,num_bts,num_users,num_ants);
kbyk_analogbeamformvecs = ones(num_chans,num_bts,num_ants);
kbyk_analogbeamformvecs_nulls = ones(num_chans,num_bts,num_ants);

chan_est_snr_db = 20;

dist_thresh = 0.5;

chan_index=1;

all_selected_cords = zeros(num_chans,num_users,2);

user_distr ='uniform'; % can be spatial or uniform

while(chan_index<=num_chans)
    
    if(user_distr == 'uniform')
        selected_user_cords = sample_points(randperm(num_choices,num_users),:);
    elseif(user_distr == 'spatial')
        selected_user_cords = zeros(num_users,2);
        for k=1:1:num_users
            y_grid_ind = floor((k-1)/4);
            x_grid_ind = mod(k-1,4);
            
            x_grid_min = (x_grid_ind)*grid_dim/4+0.5;
            x_grid_max = (x_grid_ind+1)*grid_dim/4-0.5;
            
            y_grid_min = (y_grid_ind)*grid_dim/2+0.5;
            y_grid_max = (y_grid_ind+1)*grid_dim/2-0.5;
            
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
        
%     scatter(selected_user_cords(:,1),selected_user_cords(:,2),'x')
    all_selected_cords(chan_index,:,:) = selected_user_cords;

    for j=1:1:num_bts
        curr_bts_cords = bts_cords(j,:);
        [kbyk_chans(chan_index,j,:,:), kbyk_analogbeamformvecs(chan_index,j,:,:), kbyk_analogbeamformvecs_nulls(chan_index,j,:,:)] =  chan_est_geometrical_kbyk(selected_user_cords, curr_bts_cords, chan_est_snr_db, grid_dim,  num_ants, j);
%         kbyk_chans(chan_index,j,:,:)
    end
    
    chan_index= chan_index+1;
end

legend_str= ["BTS"];

for k=1:1:num_users
    h(k+1)= scatter(all_selected_cords(:,k,1),all_selected_cords(:,k,2),'x');
    legend_str =[legend_str "User Index: "+num2str(k)];
end

h(num_users+2) = yline(grid_dim);
xline(0);
yline(0);
xline(grid_dim);
legend(h,[legend_str 'Reflector'])