%*************************************************************************************
% Negative-Binomial Fano Factor Analysis.
%*************************************************************************************
% NEURL-UA  302
% Final Project
% Author: Mark Kuang
% Email: yk2516@nyu.edu
% Creation Date: Dec 6, 2020
%*************************************************************************************

tic
close all;clear;clc;
set(0,'DefaultLineLineWidth',2);

%% 1.1 Set parameters
g = 15;
b = 0.1;
sig_tc = 5;
struevec = linspace(0,60,1000);
sprefvec = linspace(20,40,50);

f_s = NaN(length(struevec),length(sprefvec));
mean_f_s = NaN(1,length(struevec)); 

for i = 1:length(struevec)
    s = struevec(i);
    f_s(i,:) = g * exp(-(s-sprefvec).^2/2/sig_tc^2) + b;
    mean_f_s(i) = mean(f_s(i,:));
end

mean_f_s_sort = sort(mean_f_s);

figure;
plot(struevec, f_s,'.-','LineWidth',2)
set(gca,'FontSize', 18'); box off;
set(gca,'LineWidth',2); xlabel('stimuli'); ylabel('neural responses')
title('Gaussian Population Tuning Curve');

% subplot(1,2,2);plot(struevec, mean_f_s,'.-','LineWidth',2)
% set(gca,'FontSize', 18'); box off;
% set(gca,'LineWidth',2); xlabel('stimuli'); ylabel('mean responses')
% title('Mean spike count');

%% 1.2 Poisson-Gamma Encoding
num_trials = 100; 
scale_shape_s = linspace(1,81,50);

[gamma_par,modpoi] = deal(NaN(length(scale_shape_s),num_trials));
[mean_modpoi,var_modpoi]= deal(NaN(1,length(scale_shape_s)));
modpoi_data = NaN(length(scale_shape_s),num_trials,length(mean_f_s));
[mean_modpoi_data,var_modpoi_data]= deal(NaN(length(scale_shape_s),length(mean_f_s)));

for i = 1:length(scale_shape_s) % iterate over 50 different scale parameters
    
    current_iterate = i
    scale_s = scale_shape_s(i);
    %*************************************************************************************
    % test: try to see the change in shape_r
    shape_r = scale_s;
    % shape_r = 1 / scale_s;
    %*************************************************************************************
    
    for j = 1:num_trials
        % poisson-gamma stimulation without f(s)
        gamma_par(i,j) = gamrnd(shape_r,scale_s);
        modpoi(i,j) = poissrnd(gamma_par(i,j));
        
        % poisson-gamma stimulation with f(s)
        modpoi_data(i,j,:) = poissrnd(gamma_par(i,j) * mean_f_s_sort);
        
        % pars_temp = nbinfit(modpoi_data(i,j,:));

    end
        
    % poisson-gamma mean & variance with f(s)
    % gonna double check on this! this mean modpoi data & var modpoi data! 
        
    % now mean is 50x1000, for 50 increasing scale pars s, 1000 increasing
    % f(S)
    mean_modpoi_data(i,:) = mean(modpoi_data(i,:,:));
    
    % now variance is 50x1000, for 50 increasing scale pars s, 1000 increasing
    % f(S)
    var_modpoi_data(i,:) = var(modpoi_data(i,:,:));
    
end

%% plot gamma parameters

figure;
h = pcolor(gamma_par);

set(h,'Facecolor','interp')
set(h,'Linestyle','none')
set(gca,'YDir','reverse')
set(gca,'FontSize', 18);
colormap('parula')
% set(gca,'ytick',scale_shape_s)		

h = colorbar;
ylabel(h, 'mean population response bar')
xlabel('numbers of trials')
ylabel('scale parameter s indices (from 1 to 81 in 50 steps)')
title('heat map of gamma parameters')


%% 1.3 Fano factor Analysis

ff_modpoi_data = NaN(length(scale_shape_s),length(mean_f_s));
s_error = NaN(length(scale_shape_s),length(mean_f_s));

% mean_of_ff = NaN(1,length(scale_shape_s));
[beta,beta_true] = deal(NaN(1,5));

figure;

for i = 1:length(scale_shape_s)
    ff_modpoi_data(i,:) = var_modpoi_data(i,:) ./ mean_modpoi_data(i,:);
    
    s_est = (ff_modpoi_data(i,:) - 1) ./ scale_shape_s(i);
    s_error(i,:) = abs(mean_f_s_sort - s_est);
    
    if i <= 45
        subplot(1,2,1);
        plot(mean_f_s_sort,ff_modpoi_data(i,:));
        hold on;
    else
        plot(mean_f_s_sort,ff_modpoi_data(i,:));
        beta(i-45) = scale_shape_s(i);
        beta_true = string(beta);
        set(gca,'FontSize', 18'); box off;
        set(gca,'LineWidth',2); 
        xlabel('mean firing rate'); ylabel('fano factor');
        title('the growth of fano factor for growing scale parameter');
        legend(beta_true)
        hold on;
        
    end
    
end

subplot(1,2,2); hold on;
ff_end = ff_modpoi_data([46,47,48,49,50],:);

for i = 1:5
    plot(mean_f_s_sort,ff_end(i,:))
    set(gca,'FontSize', 18'); box off;
    set(gca,'LineWidth',2); 
    xlabel('mean firing rate'); ylabel('fano factor');
    title('Last five ');

end

legend(beta_true)



%% fano factor surf map
figure;
surf_fano = mesh(ff_modpoi_data);
colormap (parula (5))
c_ff = colorbar;
xlabel('Stimulus Value')
ylabel('Gamma Parameter')
zlabel('Fano Factor')
title('The Growth of Fano Factor r = s = sig_G^2')
c_ff.Label.String = 'Fano Factor Bar' ;


% %% fano factor heat map
% figure;
% h_fano = pcolor(ff_modpoi_data);
% 
% set(h_fano,'Facecolor','interp')
% set(h_fano,'Linestyle','none')
% set(gca,'YDir','reverse')
% colormap('cool')
% 
% h_fano = colorbar;
% ylabel(h_fano, 'fano factor bar')
% xlabel('stimulus values')
% ylabel('gamma parameters')
% title('the growth of fano factor with varying stimuli and gamma parameters')


%% error
figure;

for i = 2:length(scale_shape_s)
    plot(mean_f_s_sort,s_error(i,:));
    hold on
    set(gca,'FontSize', 18'); box off;
    set(gca,'LineWidth',2); 
    xlabel('mean firing rate'); ylabel('s error');
    title('error in the estimate');

end

%% error mountain
figure;
surf_error = mesh(s_error);
colormap (parula (5))
c_error = colorbar;
xlabel('Stimulus Value')
ylabel('Gamma Parameter')
zlabel('Error')
title('The Error Mountain for r = s = sig_G^2')
c_error.Label.String = 'Error Bar' ;




%% run time recording
toc
















