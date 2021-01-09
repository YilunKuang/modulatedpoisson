%*************************************************************************************


%*************************************************************************************



tic
close all;clear;clc;

%% 1 Gaussian tuning curve
%*************************************************************************************
g = 15;
%*************************************************************************************

b = 0.1;
sig_tc = 5;

%*************************************************************************************
struevec = linspace(0,60,100);
%*************************************************************************************

sprefvec = linspace(20,40,50);
f_s_nbm = NaN(length(struevec),length(sprefvec));

mean_f_s = NaN(1,length(struevec)); 

for i = 1:length(struevec)
    s = struevec(i);
    f_s_nbm(i,:) = g * exp(-(s-sprefvec).^2/2/sig_tc^2) + b;
    mean_f_s(i) = mean(f_s_nbm(i,:));
end

mean_f_s_sort = sort(mean_f_s);

figure;
plot(struevec, f_s_nbm,'.-','LineWidth',2);
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('stimuli'); ylabel('neural responses')
title('Gaussian Population Tuning Curve');

figure;
plot(struevec, mean_f_s,'.-','LineWidth',2)
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('stimuli'); ylabel('mean spike count')
title('Mean spike count of population responses for different stimuli');

%% 2 Gamma parameter set up
scale_par_s = linspace(1,81,50);
%*************************************************************************************
% shape_par_r = 1 ./ scale_par_s;
shape_par_r = scale_par_s;
%*************************************************************************************


ntrials = 100;
gam_num = NaN(length(shape_par_r),ntrials);
[mean_gam_num,var_gam_num] = deal(NaN(1,length(shape_par_r)));

% figure;
for i = 1:length(shape_par_r)
    
    for j = 1:ntrials
        gam_num(i,j) = gamrnd(shape_par_r(i),scale_par_s(i));
    end
    
%     subplot(1,2,1);
%     plot(linspace(1,10,100),gam_num(i,:),'-','LineWidth',2)
%     hold on;
%     set(gca,'FontSize', 18'); box off;
%     set(gca,'LineWidth',2); 
%     xlabel('100 trials'); 
%     ylabel('gamma values')
%     title('The fluctuations of gamma values over 100 trials');

    mean_gam_num(i) = mean(gam_num(i,:));
    var_gam_num(i) = var(gam_num(i,:));
    
    gam_num_temp = gam_num(i,:);
    gam_num_small = gam_num_temp(gam_num(i,:)<1e-5);
    gam_num_small_length = length(gam_num_small);
    
%     subplot(1,2,2);
%     ylim([-inf 1e-5])
%     plot(linspace(1,10,gam_num_small_length),gam_num_small,'-','LineWidth',2)
%     set(gca,'FontSize', 18'); box off;
%     set(gca,'LineWidth',2); 
%     hold on;
%     xlabel('100 trials'); 
%     ylabel('gamma values')
%     title('gamma values under 1e-5');

end

figure;
subplot(2,2,1); plot(linspace(1,50,50),shape_par_r,'.','LineWidth',2)
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('linspace'); ylabel('shape parameters')
title('shape parameter');

subplot(2,2,2); plot(linspace(1,50,50),scale_par_s,'.','LineWidth',2)
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('linspace'); ylabel('scale parameters')
title('scale parameter');

subplot(2,2,3); plot(linspace(1,50,50),mean_gam_num,'.','LineWidth',2)
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('shape parameter'); ylabel('mean gamma value across trials')
title('gamma mean');

subplot(2,2,4); plot(linspace(1,50,50),var_gam_num,'.','LineWidth',2)
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('shape parameter'); ylabel('variance gamma value across trials')
title('gamma variance');



%% 2.1 Heat map of gamma parameters
figure;
h_gam = pcolor(gam_num);

set(h_gam,'Facecolor','interp')
set(h_gam,'Linestyle','none')
set(gca,'YDir','reverse')
colormap('hot')

h_gam = colorbar;

ylabel(h_gam, 'mean population response bar')
xlabel('trials')
ylabel('shape r')
title('heat map of mean population response varying stimuli and gamma parameters')


%% 2.2 NBM Encoding & Decoding
% repetition set up
num_trials = 100;

% array init
[alpha,p_par,r_par,nbm_across_trial,nbm_var_across_trial] = deal(NaN(length(struevec),length(scale_par_s)));
nbm_data = NaN(length(struevec),length(scale_par_s),num_trials);
f_s_est = NaN(length(struevec),length(scale_par_s),num_trials);

for i = 1:length(struevec) % 1x100
    current_iterate_i = i
    
    % select the population response to different stimulus
    neural_population = f_s_nbm(i,:);
    mean_neural_population_sort = mean_f_s_sort(i);
    
    for j = 1:length(scale_par_s) % 1x50
        current_iterate_j = j;
        scale_s = scale_par_s(j); % 1x1
        
        % alpha = s * f(S)\Delta t        
        %*************************************************************************************
        alpha(i,j) = scale_s * mean_neural_population_sort; 
        % alpha(i,j) = scale_s * mean(neural_population); 
        %*************************************************************************************
        
        % p = 1 / (alpha+1)
        %*************************************************************************************
        p_par(i,j) = 1 / (alpha(i,j) + 1); 
        % r_par(i,j) = 1 / scale_s;
        r_par(i,j) = scale_s;
        %*************************************************************************************
            
        % repetition phase
        r_rep = repmat(r_par(i,j),[1,num_trials]); % 1 x num_trials
        p_rep = repmat(p_par(i,j),[1,num_trials]); % 1 x num_trials
        
        % mean spi
        nbm_data(i,j,:) = nbinrnd(r_rep,p_rep); % 1 x num_trials
        
        nbm_across_trial(i,j) = mean(nbm_data(i,j,:));
        nbm_var_across_trial(i,j) = var(nbm_data(i,j,:));
        
        for k = 1:num_trials % 1x100
            current_iterate_k = k;
            
            % MLE for NBM parameters
            pars_temp = nbinfit(nbm_data(i,j,:));
            % f(S) decoding estimate
            f_s_est(i,j,k) = (pars_temp(1) * (1 - pars_temp(2))) / pars_temp(2);

%             mean_temp = mean(nbm_data(:,k));
%             var_temp = var(nbm_data(:,k));
%             
%             ff_temp = var_temp / mean_temp;
%             
%             if ff_temp <= 1
%                 f_s_est(i,j,k) = 0;
%             else
%                 % MLE for NBM parameters
%                 pars_temp = nbinfit(nbm_data(:,k));
%                 % f(S) decoding estimate
%                 f_s_est(i,j,k) = (pars_temp(1) * (1 - pars_temp(2))) / pars_temp(2);
%                 
%             end
            
            
        end

    end
    
    % nan-check -> seems like there is no nan in here. Then what's the
    % problem?
    % sum(nan_check(:))
    
end

%% 

% heat map of mean population responses
figure;
h_nbm = pcolor(nbm_across_trial);
set(h_nbm,'Facecolor','interp')
set(h_nbm,'Linestyle','none')
set(gca,'YDir','reverse')
colormap('hot')

h_nbm = colorbar;
ylabel(h_nbm, 'mean population response bar')
xlabel('gamma parameters')
ylabel('stimulus values')
title('heat map of mean population response varying stimuli and gamma parameters')

% % population response against gamma 
% nbm_across_stim = mean(nbm_across_trial);
% xx = linspace(1,50,50);
% figure;
% plot(xx, nbm_across_stim,'.-','LineWidth',2);
% set(gca,'FontSize', 18'); box off;
% set(gca,'LineWidth',2); xlabel('shape parameter r'); ylabel('population response averaged over stimuli ')
% title('the growth of population response against gamma');

%% fano factor analysis
% ff_modpoi_data = NaN(length(scale_par_s),length(struevec));
% 
% nbm_var_across_trial_T = transpose(nbm_var_across_trial);
% nbm_across_trial_T = transpose(nbm_across_trial);
% 
% figure;
% for i = 1:length(scale_par_s)
%     
%     ff_modpoi_data(i,:) = nbm_var_across_trial_T(i,:) ./ nbm_across_trial_T(i,:);
%     % s_est = (ff_modpoi_data(i,:) - 1) ./ scale_par_s(i);
%     
%     if i < 40
%         plot(mean_f_s_sort,ff_modpoi_data(i,:));
%         hold on
%         set(gca,'FontSize', 18'); box off;
%         set(gca,'LineWidth',2); 
%         xlabel('mean firing rate'); ylabel('fano factor');
%         title('the growth of fano factor for growing scale parameter');
%     else
%         beta = string(scale_par_s(i));
%         plot(mean_f_s_sort,ff_modpoi_data(i,:));
%         legend('beta='+beta)
%         hold on
%         set(gca,'FontSize', 18'); box off;
%         set(gca,'LineWidth',2); 
%         xlabel('mean firing rate'); ylabel('fano factor');
%         title('the growth of fano factor for growing scale parameter');
% 
%     end
%     
% end



%% decoding
tic
% load('f_s_est.mat')

[nan_count,mean_f_s_est,bias_ML,var_ML] = deal(NaN(length(struevec),length(scale_par_s)));

for i = 1:length(struevec) % 1x100
    s_compared = mean_f_s(i);
    
    for j = 1:length(scale_par_s) % 1x50
        nan_count(i,j) = sum(isnan(f_s_est(i,j,:)));
        
        if nan_count(i,j) == 0
            mean_f_s_est(i,j) = mean(f_s_est(i,j,:));
            bias_ML(i,j) = mean_f_s_est(i,j) - s_compared;
            var_ML(i,j) = var(f_s_est(i,j,:));
        else
            mean_f_s_est(i,j) = NaN;
            bias_ML(i,j) = NaN;
            var_ML(i,j) = NaN;
        end
        
                
    end
    
end

toc

%% KNN impute

bias_ML_KNN = knnimpute(bias_ML);
var_ML_KNN = knnimpute(var_ML);

%% bias heatmap

% heat map of mean population responses
figure;
h_ML_bias = pcolor(bias_ML_KNN);
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); 
set(h_ML_bias,'Facecolor','interp')
set(h_ML_bias,'Linestyle','none')
set(gca,'YDir','reverse')
colormap('parula')

h_ML_bias = colorbar;
ylabel(h_ML_bias, 'decoding bias bar')
xlabel('gamma parameter (from 1 to 81 in steps of 50)')
ylabel('f(S) (from 0.1 to 8.9 in steps of 100)')
title('decoding bias for r = s = sig_G^2')

%% variance heatmap

% heat map of mean population responses
figure;
h_ML_var = pcolor(var_ML_KNN);
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); 
set(h_ML_var,'Facecolor','interp')
set(h_ML_var,'Linestyle','none')
set(gca,'YDir','reverse')
colormap('parula')

h_ML_var = colorbar;
ylabel(h_ML_var, 'decoding variance bar')
xlabel('gamma parameter (from 1 to 81 in steps of 50)')
ylabel('stimulus values')
title('decoding variance')

%% Two-sample Kolmogorov-Smirnov test
% test_val = kruskalwallis(bias_ML_KNN);

[h_ks,p_ks] = deal(NaN(length(shape_par_r),length(shape_par_r)));

for i = 1:length(shape_par_r) 
    gam_sample_i = bias_ML_KNN(:,i);
    
    for j = 1:length(shape_par_r)
        gam_sample_j = bias_ML_KNN(:,j);
        [h_ks(i,j),p_ks(i,j)] = kstest2(gam_sample_i,gam_sample_j);
        
    end
    
end

%% Plotting

p_check = sum(p_ks < 0.05);


figure;
h_KS = pcolor(p_ks);
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); 
set(h_KS,'Facecolor','interp')
set(h_KS,'Linestyle','none')
set(gca,'YDir','reverse')
colormap('parula')

h_KS = colorbar;
ylabel(h_KS, 'p-value bar')
xlabel('sample i for i = 1:50')
ylabel('sample j for j = 1:50')
title('p-value map for r = s = sig_G^2')


    
    


%% Runtime recording
toc




