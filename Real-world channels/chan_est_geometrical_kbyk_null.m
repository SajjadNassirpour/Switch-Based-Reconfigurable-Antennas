% Example usage, user_locs = [[3,3],[4,4],[-1,2]], bas_station_loc = [0,0],
% chan_est_snr_db = 40
% chan_est_geometrical([[3,3];[4,4];[-1,2]],[0,0],40);

function [chan_est,chan_est_nonoise,analogbeamforming_vec,analogbeamforming_vec_null,analogbeamforming_vec_null_multipath] = chan_est_geometrical_kbyk_null(user_locs, base_station_loc, chan_est_snr_db, grid_dim,  num_ants, support_user_index, enable_multipath, freq, rx_ant_sep_lambda, multipath_loss)
  % Default params 
  if nargin <= 6
    enable_multipath = true;
    freq = 2.4e9;
    rx_ant_sep_lambda = 0.25;
    multipath_loss = 1;
  end
  
  reflector_slope_intercepts = [[0,grid_dim];[90,grid_dim];[90,0];[0,0]];
  [num_reflectors,~] = size(reflector_slope_intercepts); 
  
  wav = 3e8/freq;
  [num_users,~] = size(user_locs);
  rx_ant_sep = rx_ant_sep_lambda*wav;
  
  bts_cords = get_rx_ant_cords(base_station_loc,1,num_ants, rx_ant_sep, 0);
  
  chan_est = zeros([num_users,num_ants]);
  chan_est_nonoise = zeros([num_users,num_ants]);
  analogbeamforming_vec_all = zeros([num_users,num_ants]);
  analogbeamforming_vec_all_withmp = zeros([num_users,num_ants]);
  
  analogbeamforming_vec_null = zeros([num_ants]);
  analogbeamforming_vec_null_multipath = zeros([num_ants]);
  

  for j = 1:1:num_users
      tiled_user_locs = repmat(user_locs(j,:),[num_ants,1]);
      curr_bts_cords = squeeze(bts_cords(1,:,:));
      % Norm of order 2, along dimension 2
      tx_rx_dists = vecnorm(tiled_user_locs-curr_bts_cords,2,2).';
      avg_dir_loss = 1/mean(tx_rx_dists);
      if(j==1)
          norm_fac = avg_dir_loss;
      end
      chan_mat = exp(-1j*(2*pi*freq/3e8)*(tx_rx_dists))*avg_dir_loss;
      dir_path_chan = chan_mat;
      
      
      if(enable_multipath)
          for reflector_index=1:1:num_reflectors
            [reflected_distances, path_loss_facs] = get_reflected_dists(tiled_user_locs,curr_bts_cords,reflector_slope_intercepts(reflector_index,:));
            reflected_distances = reflected_distances.';
            avg_refl_loss = 1/(mean(path_loss_facs));
            chan_mat = chan_mat+exp(-1j*(2*pi*freq/3e8)*(reflected_distances))*avg_refl_loss*multipath_loss;
          end
      end
      
      % Noise added in reference to user 1
      chan_mat_txpow = chan_mat/norm_fac;
      dirpath_chan_mat_txpow = dir_path_chan/norm_fac;
      
      if(j==1)
        noise_pow = 10^((-chan_est_snr_db)/10);
        noise_std_dev = sqrt(noise_pow);
      end

      noise_sig = normrnd(0,noise_std_dev/sqrt(2),1,num_ants)+1j*normrnd(0,noise_std_dev/sqrt(2),1,num_ants);
      
      chan_est_nonoise(j,:) = chan_mat_txpow;
      chan_est(j,:) = chan_mat_txpow+noise_sig;
      analogbeamforming_vec_all(j,:) = dirpath_chan_mat_txpow+noise_sig;
      analogbeamforming_vec_all_withmp(j,:) = chan_mat_txpow+noise_sig;

  end
   
  support_vec = squeeze(analogbeamforming_vec_all(support_user_index,:));
  
  interferers = 1:1:num_users;
  interferers(support_user_index)=[];
  
  intrf_mat = analogbeamforming_vec_all(interferers,:);
  
  null_mat = null(intrf_mat);
  
  null_comb_vec = conj(support_vec*null_mat);
  
  analogbeamforming_vec_null = null_comb_vec*(null_mat.');
  analogbeamforming_vec_null = analogbeamforming_vec_null/norm(analogbeamforming_vec_null);
  
  analogbeamforming_vec = conj(support_vec);
  analogbeamforming_vec = analogbeamforming_vec/norm( analogbeamforming_vec);
  
  support_vec_mp = squeeze(analogbeamforming_vec_all_withmp(support_user_index,:));
  intrf_mat_mp = analogbeamforming_vec_all_withmp(interferers,:);
  
  null_mat = null(intrf_mat_mp);
  
  
  null_comb_vec = conj(support_vec_mp*null_mat);
  
  analogbeamforming_vec_null_multipath = null_comb_vec*(null_mat.');
  analogbeamforming_vec_null_multipath = analogbeamforming_vec_null_multipath/norm(analogbeamforming_vec_null_multipath);
  
end

function ant_cords = get_rx_ant_cords(base_station_loc, num_rfc,num_ants, rx_ant_sep, rfc_step_fac)
    num_ants_onecol = 2;
    num_ants_onerow = floor(num_ants/num_ants_onecol);
    
    ant_cords = zeros([num_rfc,num_ants,2]);
    
    dist_vec_row = 0:rx_ant_sep:(num_ants_onerow-1)*rx_ant_sep;
    dist_vec_col = 0:rx_ant_sep:(num_ants_onecol-1)*rx_ant_sep;
    rf_chan_sep = num_ants_onerow*rx_ant_sep*rfc_step_fac;
    additional_rfc_incr = rf_chan_sep*[1,0];
    
    for i=1:1:num_rfc
        [X,Y] = meshgrid(dist_vec_row,dist_vec_col);
        cart_prod = [X(:) Y(:)];
        ant_cords(i,:,:) = base_station_loc+(i-1)*additional_rfc_incr+ cart_prod;
    end
    
end

function [reflected_dists, path_loss_facs] = get_reflected_dists(tiled_user_locs,dest_pts,reflector_slope_intercept)
    image_pts = zeros(size(dest_pts));
    int_pts = zeros(size(dest_pts));
    if(reflector_slope_intercept(1)~=90)
        m = tan(reflector_slope_intercept(1)*pi/180);
        c = reflector_slope_intercept(2);
        
        x1 = tiled_user_locs(1,1);
        y1 = tiled_user_locs(1,2);
        
        image_pts(:,1) = -2*m*(m*x1-y1+c)/(m^2+1)+x1;
        image_pts(:,2) = 2*(m*x1-y1+c)/(m^2+1)+y1;
        reflected_dists = vecnorm(image_pts-dest_pts,2,2);
        
        m1 = (image_pts(:,2)-dest_pts(:,2))./(image_pts(:,1)-dest_pts(:,1));
        c1 = image_pts(:,2)-m1.*image_pts(:,1);
        
        
        int_pts(:,1) = (c1-c)./(m-m1);
        int_pts(:,2) = m*int_pts(:,1)+c;
        
        path_loss_facs=vecnorm(int_pts-image_pts,2,2).*vecnorm(int_pts-dest_pts,2,2);
        
    else

        c = reflector_slope_intercept(2);
        
        x1 = tiled_user_locs(1,1);
        y1 = tiled_user_locs(1,2);
        
        image_pts(:,1) = 2*c-x1;
        image_pts(:,2) = y1;
        
        reflected_dists = vecnorm(image_pts-dest_pts,2,2);

        m1 = (image_pts(:,2)-dest_pts(:,2))./(image_pts(:,1)-dest_pts(:,1));
        c1 = image_pts(:,2)-m1.*image_pts(:,1);
        
        int_pts(:,1) = c;
        int_pts(:,2) = m1.*int_pts(:,1)+c1;

        
        path_loss_facs=vecnorm(int_pts-image_pts,2,2).*vecnorm(int_pts-dest_pts,2,2);
    end
end




