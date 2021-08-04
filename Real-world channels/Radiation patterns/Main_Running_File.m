% ***** Radiation Patterns ***** %
clc
clear 
close all

rng(15);

disable_switched_app = false;
K=4;  % Number of users (BTSs)
num_users=K;
M=64; % NUmber of elements in array antenna

est_SNR_mat=20;  % SNR(dB) in channel estimation process

iters=1;  % Number of iterations

pow_mat=10;  % SNR interval

% Optimization parameters in our approach
mu=0.0001;
r=10^5;
divide_rate=2;

% Sum-rate with Analog beamforming
sum_rate_analog=zeros(length(est_SNR_mat),1,length(pow_mat));  

% Sum-rate with TDMA+Analog beamforming
sum_rate_tdma_analog=zeros(length(est_SNR_mat),1,length(pow_mat));

% Sum-rate with Beamnulling+Analog beamforming
sum_rate_beam_null=zeros(length(est_SNR_mat),1,length(pow_mat));

% Sum-rate with Beamnulling(multipath)+Analog beamforming
sum_rate_beam_null_multipath=zeros(length(est_SNR_mat),1,length(pow_mat));

% Sum-rate with MMSE
sum_rate_mmse=zeros(length(est_SNR_mat),2,length(pow_mat));

% Sum-rate with our approach
sum_rate_our=zeros(length(est_SNR_mat),1,length(pow_mat));


for jj=1:length(est_SNR_mat)
    
    est_SNR=est_SNR_mat(jj);
    
    % Generate the channels
    [H_noisy,H_nonoise,beam_vector, null_vector, null_multipath_vector,bts_cords,all_selected_cords,grid_dim,bts_cord_ants,refl_cords]=real_world_channel_beam_null(K,M,iters,est_SNR, 'uniform');
 
    
    % MMSE
    sum_rate_mmse(jj,:,:)=MMSE_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise);
    mmse_iter=squeeze(sum_rate_mmse(jj,:,:))
    
    % Analog Beamforming
    sum_rate_analog(jj,:,:)=Analog_app_complex(K,M,iters,pow_mat,H_nonoise,beam_vector);
    analog_iter=squeeze(sum_rate_analog(jj,:,:))'
    
    % TDMA + Analog Beamforming
    sum_rate_tdma_analog(jj,:,:)=TDMA_Analog_app_complex(K,M,iters,pow_mat,H_nonoise,beam_vector);
    tdma_analog_iter=squeeze(sum_rate_tdma_analog(jj,:,:))'
    
    % Beamnulling + Analog Beamforming
    sum_rate_beam_null(jj,:,:)=Beam_null_app_complex(K,M,iters,pow_mat,H_nonoise,null_vector);
    beam_null_iter=squeeze(sum_rate_beam_null(jj,:,:))'

    % noise_cond
    % false: Noiseless channels
    % true: Noisy channels 
    noise_cond=true; 
    
    % Switch_mat : Optimal switch configuration 
    if(~disable_switched_app)
        [sum_rate_our(jj,:,:), Switch_mat]=Binary_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise,mu,r,divide_rate,noise_cond);
        our_iter=squeeze(sum_rate_our(jj,:,:))'
        switch_vec = Switch_mat;
    end
    
    % ---- Radiation Patterns ---- %
    freq=2.4e9;
    beam_test = squeeze(beam_vector);
    figure(1)
    subplot(2,2,1)
    
    % -- Analog Beamforming -- %
    leg=[];
    leg = [leg, scatter(bts_cords(:,1),bts_cords(:,2),100,'s','filled')];
    hold on

    ylim([-2,grid_dim+2]);
    xlim([-2,grid_dim+2]);

    legend_str= ["BTS"];

    for k=1:1:num_users
        leg(k+1)= scatter(all_selected_cords(:,k,1),all_selected_cords(:,k,2),100,'x','LineWidth',3);
        legend_str =[legend_str "User : "+num2str(k)];
    end

    leg(num_users+2) = yline(grid_dim);
    xline(0);
    yline(0);
    xline(grid_dim);
    legend_str =[legend_str 'Reflector'];
    first_bool = true;
    for bts_ind=1:1:K
        beam_vec_bts = squeeze(beam_test(bts_ind,:)).';
        curr_cord = squeeze(bts_cords(bts_ind,:));
        curr_user = squeeze(all_selected_cords(1,bts_ind,:));
        curr_ant_cords = squeeze(bts_cord_ants(1,bts_ind,:,:));
        
        test_angs = 0:1:359;
        test_cords = zeros([numel(test_angs),2]);
        test_radii = grid_dim/4;
        test_cords(:,1) = test_radii*cos(test_angs.*pi/180);
        test_cords(:,2) = test_radii*sin(test_angs.*pi/180);
        rad_pattern = zeros([1,numel(test_angs)]);
        
        for test_ind = 1:1:numel(test_angs)
            curr_test_cord = squeeze(test_cords(test_ind,:))+curr_cord;
            distances = vecnorm(repmat(curr_test_cord,[M,1])-curr_ant_cords,2,2);
            phase_vec = exp(-1j*((2*pi/3e8)*freq)*distances);
            rad_pattern(test_ind) = abs(sum(phase_vec.*beam_vec_bts))^2;
        end
        
        rad_pattern = rad_pattern./max(rad_pattern);
        test_cords = zeros([numel(test_angs),2]);
        test_radii = grid_dim/6;
        test_cords(:,1) = test_radii*cos(test_angs.*pi/180).*rad_pattern;
        test_cords(:,2) = test_radii*sin(test_angs.*pi/180).*rad_pattern;
        plot_cords = repmat(curr_cord,[numel(test_angs),1])+test_cords;
        
        
        if(bts_ind==1)
            leg(num_users+3)=plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2);
            leg(num_users+4)=line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5);
            legend_str =[legend_str ,'Beam Pattern','User'];
        else
            plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2)
            line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5)
        end
        
        interferers = 1:1:num_users;
        interferers(bts_ind)=[]; 
    end
    grid on
    title('$\textbf{Analog Beamforming}$','Interpreter','latex')
    xlabel('x','Interpreter','latex')
    ylabel('y','Interpreter','latex')
    grid on
    set(gcf,'color','w');
    set(gca,'TickLabelInterpreter','latex')
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'Interpreter','latex','fontSize',15)
    set(gca,'fontsize',15)

    % -- Beamnulling + Analog Beamforming -- %
    beam_test = squeeze(null_vector);
    figure(1)
    subplot(2,2,2)
    leg=[];

    leg = [leg, scatter(bts_cords(:,1),bts_cords(:,2),100,'s','filled')];
    hold on

    ylim([-2,grid_dim+2]);
    xlim([-2,grid_dim+2]);

    legend_str= ["BTS"];

    for k=1:1:num_users
        leg(k+1)= scatter(all_selected_cords(:,k,1),all_selected_cords(:,k,2),100,'x','LineWidth',3);
        legend_str =[legend_str "User : "+num2str(k)];
    end

    leg(num_users+2) = yline(grid_dim);
    xline(0);
    yline(0);
    xline(grid_dim);
    legend_str =[legend_str 'Reflector'];
    first_bool = true;
    for bts_ind=1:1:K
        beam_vec_bts = squeeze(beam_test(bts_ind,:)).';
        curr_cord = squeeze(bts_cords(bts_ind,:));
        curr_user = squeeze(all_selected_cords(1,bts_ind,:));
        curr_ant_cords = squeeze(bts_cord_ants(1,bts_ind,:,:));
        
        test_angs = 0:1:359;
        test_cords = zeros([numel(test_angs),2]);
        test_radii = grid_dim/4;
        test_cords(:,1) = test_radii*cos(test_angs.*pi/180);
        test_cords(:,2) = test_radii*sin(test_angs.*pi/180);
        rad_pattern = zeros([1,numel(test_angs)]);
        
        for test_ind = 1:1:numel(test_angs)
            curr_test_cord = squeeze(test_cords(test_ind,:))+curr_cord;
            distances = vecnorm(repmat(curr_test_cord,[M,1])-curr_ant_cords,2,2);
            phase_vec = exp(-1j*((2*pi/3e8)*freq)*distances);
            rad_pattern(test_ind) = abs(sum(phase_vec.*beam_vec_bts))^2;
        end
        
        rad_pattern = rad_pattern./max(rad_pattern);
        test_cords = zeros([numel(test_angs),2]);
        test_radii = grid_dim/6;
        test_cords(:,1) = test_radii*cos(test_angs.*pi/180).*rad_pattern;
        test_cords(:,2) = test_radii*sin(test_angs.*pi/180).*rad_pattern;
        plot_cords = repmat(curr_cord,[numel(test_angs),1])+test_cords;
        
        
        if(bts_ind==1)
            leg(num_users+3)=plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2);
            leg(num_users+4)=line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5);
            legend_str =[legend_str ,'Beam Pattern','User'];
        else
            plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2)
            line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5)
        end
        
        interferers = 1:1:num_users;
        interferers(bts_ind)=[];
        
        
        for intf_idx=interferers
            intf_user = squeeze(all_selected_cords(1,intf_idx,:));
            if(first_bool)
                leg(num_users+5)=line([curr_cord(1),intf_user(1)],[curr_cord(2),intf_user(2)],'Color','red','LineStyle',':','LineWidth',0.5);
                legend_str = [legend_str 'Interferers'];
                first_bool = false;
            else
                line([curr_cord(1),intf_user(1)],[curr_cord(2),intf_user(2)],'Color','red','LineStyle',':','LineWidth',0.5)
            end
        end
    end
    grid on
    title('$\textbf{Analog Beamforming+Beamnulling}$','Interpreter','latex')
    xlabel('x','Interpreter','latex')
    ylabel('y','Interpreter','latex')
    grid on
    set(gcf,'color','w');
    set(gca,'TickLabelInterpreter','latex')
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'Interpreter','latex','fontSize',15)
    set(gca,'fontsize',15)
    
    % -- MMSE -- %
    beam_test = squeeze(null_multipath_vector);
    figure(1)
    subplot(2,2,3)
    leg=[];
    leg = [leg, scatter(bts_cords(:,1),bts_cords(:,2),100,'s','filled')];
    hold on
    ylim([-2,grid_dim+2]);
    xlim([-2,grid_dim+2]);
    legend_str= ["BTS"];

    for k=1:1:num_users
        leg(k+1)= scatter(all_selected_cords(:,k,1),all_selected_cords(:,k,2),100,'x','LineWidth',3);
        legend_str =[legend_str "User : "+num2str(k)];
    end

    leg(num_users+2) = yline(grid_dim);
    xline(0);
    yline(0);
    xline(grid_dim);
    legend_str =[legend_str 'Reflector'];
    first_bool = true;
    first_bool1 = true;
    first_bool2 = true;
    for bts_ind=1:1:K
        beam_vec_bts = squeeze(beam_test(bts_ind,:)).';
        curr_cord = squeeze(bts_cords(bts_ind,:));
        curr_user = squeeze(all_selected_cords(1,bts_ind,:));
        curr_ant_cords = squeeze(bts_cord_ants(1,bts_ind,:,:));
        
        test_angs = 0:1:359;
        test_cords = zeros([numel(test_angs),2]);
        test_radii = grid_dim/4;
        test_cords(:,1) = test_radii*cos(test_angs.*pi/180);
        test_cords(:,2) = test_radii*sin(test_angs.*pi/180);
        rad_pattern = zeros([1,numel(test_angs)]);
        
        for test_ind = 1:1:numel(test_angs)
            curr_test_cord = squeeze(test_cords(test_ind,:))+curr_cord;
            distances = vecnorm(repmat(curr_test_cord,[M,1])-curr_ant_cords,2,2);
            phase_vec = exp(-1j*((2*pi/3e8)*freq)*distances);
            rad_pattern(test_ind) = abs(sum(phase_vec.*beam_vec_bts))^2;
        end
        
        rad_pattern = rad_pattern./max(rad_pattern);
        test_cords = zeros([numel(test_angs),2]);
        test_radii = grid_dim/6;
        test_cords(:,1) = test_radii*cos(test_angs.*pi/180).*rad_pattern;
        test_cords(:,2) = test_radii*sin(test_angs.*pi/180).*rad_pattern;
        plot_cords = repmat(curr_cord,[numel(test_angs),1])+test_cords;
        
        
        if(bts_ind==1)
            leg(num_users+3)=plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2);
            leg(num_users+4)=line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5);
            legend_str =[legend_str ,'Beam Pattern','User'];
        else
            plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2)
            line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5)
        end
        
        interferers = 1:1:num_users;
        interferers(bts_ind)=[];
        
        
        for intf_idx=interferers
            intf_user = squeeze(all_selected_cords(1,intf_idx,:));
            if(first_bool)
                leg(num_users+5)=line([curr_cord(1),intf_user(1)],[curr_cord(2),intf_user(2)],'Color','red','LineStyle',':','LineWidth',0.5);
                legend_str = [legend_str 'Interferers'];
                first_bool = false;
            else
                line([curr_cord(1),intf_user(1)],[curr_cord(2),intf_user(2)],'Color','red','LineStyle',':','LineWidth',0.5)
            end
        end
        
        refl_on_off = [[0,0,1,0];[0,1,0,0];[0,0,0,0];[1,0,0,0]];
        for usr_ind=1:1:K
            intf_user = squeeze(all_selected_cords(1,usr_ind,:));
            for refl_ind=1:1:4
                if(usr_ind==bts_ind && refl_on_off(usr_ind,refl_ind)==1)
                    col1='cyan';
                    lwid=1.5;
                    if(first_bool1)
                        leg(num_users+6)=line([curr_cord(1),refl_cords(1,bts_ind,usr_ind,refl_ind,1)],[curr_cord(2),refl_cords(1,bts_ind,usr_ind,refl_ind,2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                        legend_str = [legend_str 'Multipath Reflections'];
                        line([refl_cords(1,bts_ind,usr_ind,refl_ind,1),intf_user(1)],[refl_cords(1,bts_ind,usr_ind,refl_ind,2),intf_user(2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                        first_bool1=false;
                    else
                        line([curr_cord(1),refl_cords(1,bts_ind,usr_ind,refl_ind,1)],[curr_cord(2),refl_cords(1,bts_ind,usr_ind,refl_ind,2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                        line([refl_cords(1,bts_ind,usr_ind,refl_ind,1),intf_user(1)],[refl_cords(1,bts_ind,usr_ind,refl_ind,2),intf_user(2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                    end
                end
            end
        end
    end
    grid on
    title('$\textbf{MMSE}$','Interpreter','latex')
    xlabel('x','Interpreter','latex')
    ylabel('y','Interpreter','latex')
    grid on
    set(gcf,'color','w');
    set(gca,'TickLabelInterpreter','latex')
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'Interpreter','latex','fontSize',15)
    set(gca,'fontsize',15)
    
    % -- Our approach -- %
    if(~disable_switched_app)
        beam_test = squeeze(switch_vec);
        figure(1)
        subplot(2,2,4)
        leg=[];
        leg = [leg, scatter(bts_cords(:,1),bts_cords(:,2),100,'s','filled')];
        hold on
        ylim([-2,grid_dim+2]);
        xlim([-2,grid_dim+2]);
        legend_str= ["BTS"];
        for k=1:1:num_users
            leg(k+1)= scatter(all_selected_cords(:,k,1),all_selected_cords(:,k,2),100,'x','LineWidth',3);
            legend_str =[legend_str "User "+num2str(k)];
        end
        leg(num_users+2) = yline(grid_dim);
        xline(0);
        yline(0);
        xline(grid_dim);
        legend_str =[legend_str 'Reflector'];
        first_bool = true;
        first_bool1 = true;
        first_bool2 = true;
        for bts_ind=1:1:K
            beam_vec_bts = squeeze(beam_test(bts_ind,:)).';
            curr_cord = squeeze(bts_cords(bts_ind,:));
            curr_user = squeeze(all_selected_cords(1,bts_ind,:));
            curr_ant_cords = squeeze(bts_cord_ants(1,bts_ind,:,:));

            test_angs = 0:1:359;
            test_cords = zeros([numel(test_angs),2]);
            test_radii = grid_dim/4;
            test_cords(:,1) = test_radii*cos(test_angs.*pi/180);
            test_cords(:,2) = test_radii*sin(test_angs.*pi/180);
            rad_pattern = zeros([1,numel(test_angs)]);

            for test_ind = 1:1:numel(test_angs)
                curr_test_cord = squeeze(test_cords(test_ind,:))+curr_cord;
                distances = vecnorm(repmat(curr_test_cord,[M,1])-curr_ant_cords,2,2);
                phase_vec = exp(-1j*((2*pi/3e8)*freq)*distances);
                rad_pattern(test_ind) = abs(sum(phase_vec.*beam_vec_bts))^2;
            end

            rad_pattern = rad_pattern./max(rad_pattern);
            test_cords = zeros([numel(test_angs),2]);
            test_radii = grid_dim/6;
            test_cords(:,1) = test_radii*cos(test_angs.*pi/180).*rad_pattern;
            test_cords(:,2) = test_radii*sin(test_angs.*pi/180).*rad_pattern;
            plot_cords = repmat(curr_cord,[numel(test_angs),1])+test_cords;


            if(bts_ind==1)
                leg(num_users+3)=plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2);
                leg(num_users+4)=line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5);
                legend_str =[legend_str ,'Beam Pattern','Direct Signal'];
            else
                plot(plot_cords(:,1),plot_cords(:,2),'Color','black','LineStyle','-','LineWidth',2)
                line([curr_cord(1),curr_user(1)],[curr_cord(2),curr_user(2)],'Color','blue','LineStyle','--','LineWidth',1.5)
            end

            interferers = 1:1:num_users;
            interferers(bts_ind)=[];


            for intf_idx=interferers
                intf_user = squeeze(all_selected_cords(1,intf_idx,:));
                if(first_bool)
                    leg(num_users+5)=line([curr_cord(1),intf_user(1)],[curr_cord(2),intf_user(2)],'Color','red','LineStyle',':','LineWidth',0.5);
                    legend_str = [legend_str 'Interferers'];
                    first_bool = false;
                else
                    line([curr_cord(1),intf_user(1)],[curr_cord(2),intf_user(2)],'Color','red','LineStyle',':','LineWidth',0.5)
                end
            end
            refl_on_off = [[1,0,0,0];[1,0,0,0];[0,0,0,0];[0,0,1,0]];
            for usr_ind=1:1:K
                intf_user = squeeze(all_selected_cords(1,usr_ind,:));
                for refl_ind=1:1:4
                    if(usr_ind==bts_ind && refl_on_off(usr_ind,refl_ind)==1)
                        col1='cyan';
                        lwid=1.5;
                        if(first_bool1)
                            leg(num_users+6)=line([curr_cord(1),refl_cords(1,bts_ind,usr_ind,refl_ind,1)],[curr_cord(2),refl_cords(1,bts_ind,usr_ind,refl_ind,2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                            legend_str = [legend_str 'Multipath Signal'];
                            line([refl_cords(1,bts_ind,usr_ind,refl_ind,1),intf_user(1)],[refl_cords(1,bts_ind,usr_ind,refl_ind,2),intf_user(2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                            first_bool1=false;
                        else
                            line([curr_cord(1),refl_cords(1,bts_ind,usr_ind,refl_ind,1)],[curr_cord(2),refl_cords(1,bts_ind,usr_ind,refl_ind,2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                            line([refl_cords(1,bts_ind,usr_ind,refl_ind,1),intf_user(1)],[refl_cords(1,bts_ind,usr_ind,refl_ind,2),intf_user(2)],'Color',col1,'LineStyle',':','LineWidth',lwid);
                        end
                    end
                end
            end

        end
        grid on
        title('$\textbf{Our proposed approach}$','Interpreter','latex')
        xlabel('x','Interpreter','latex')
        ylabel('y','Interpreter','latex')
        
        grid on
        set(gcf,'color','w');
        set(gca,'TickLabelInterpreter','latex')
        figureHandle = gcf;
        set(findall(figureHandle,'type','text'),'Interpreter','latex','fontSize',15)
        set(gca,'fontsize',15)
    end

    H9 = legend(leg,legend_str,'Interpreter','latex','fontSize',15);
    set(H9,'Interpreter','latex','fontSize',15);
end



