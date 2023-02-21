%% General Script for the 10000 MC simulations
% 90% uncertainty bounds
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

%%
clear
load everything_m05_ihacres.mat
model = model_list{1};
disp(1)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m07_gr4j.mat
model = model_list{1};
disp(2)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m14_topmodel.mat
model = model_list{1};
disp(3)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m26_flexi.mat
model = model_list{1};
disp(4)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m35_mopex5.mat
model = model_list{1};
disp(5)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m36_modhydrolog.mat
model = model_list{1};
disp(6)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m37_hbv.mat
model = model_list{1};
disp(7)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end

clear
load everything_m44_echo.mat
model = model_list{1};
disp(8)
fifth_percentile = ((5/100)*numSample)+1;
nintyfifth_percentile = ((95/100)*numSample)-1;

for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    for j = 1:length(t)
        Daily_timestep_mean{j,1} = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    Daily_timestep_mean = cell2mat(Daily_timestep_mean);

    % Samples' Performance criteria
    for j = 1:numSample
        qsim = Q(j,:);
        LL(1,j)       = of_loglikelihood(qobs,qsim,warmup);
        aic(1,j)      = of_aic(qobs,qsim,warmup,numPar);
        bic(1,j)      = of_bic(qobs,qsim,warmup,numPar);
        KGE(1,j)      = of_KGE(qobs,qsim,warmup);
        KGEi(1,j)     = of_inverse_KGE(qobs,qsim,warmup);
        KGEm(1,j)     = of_mean_hilo_KGE(qobs,qsim,warmup);
        NSEf(1,j)     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSEi(1,j)     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        NSElog(1,j)   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
        RMSEE(1,j)     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)    = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = cell2mat(input_parameter(1,1));
    [~,~,w_listLL]                 = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]          = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]          = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]          = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]        = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]        = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]          = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]        = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]        = efficiency_rank(RMSE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]      = efficiency_rank(NRMSE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]              = efficiency_rank(R(1,:),numSample,'R',Par);

    % Ensemble Mean Performance Criteria
    qsim = Daily_timestep_mean.';
    simLL       = of_loglikelihood(qobs,qsim,warmup);
    simaic      = of_aic(qobs,qsim,warmup,numPar);
    simbic      = of_bic(qobs,qsim,warmup,numPar);
    simKGE      = of_KGE(qobs,qsim,warmup);
    simKGEi     = of_inverse_KGE(qobs,qsim,warmup);
    simKGEm     = of_mean_hilo_KGE(qobs,qsim,warmup);
    simNSEf     = NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSEi     = inverse_NSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNSElog   = logNSE(qobs(warmup:end,1),qsim(1,warmup:end).');
    simRMSEE     = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE    = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
    simR        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');

    b{1,1} = 'KGE';
    b{1,2} = 'KGEi';
    b{1,3} = 'KGEm';
    b{1,4} = 'NSE';
    b{1,5} = 'NSEi';
    b{1,6} = 'NSElog';
    b{1,7} = 'RMSE';
    b{1,8} = 'NRMSE';
    b{1,9} = 'R';
    b{1,10} = 'LL';
    b{1,11} = 'AIC';
    b{1,12} = 'BIC';
    b{2,1} = simKGE;
    b{2,2} = simKGEi;
    b{2,3} = simKGEm;
    b{2,4} = simNSEf;
    b{2,5} = simNSEi;
    b{2,6} = simNSElog;
    b{2,7} = simRMSEE;
    b{2,8} = simNRMSEE;
    b{2,9} = simR;
    b{2,10} = simLL;
    b{2,11} = simaic;
    b{2,12} = simbic;

    writecell(b,'Ensemble_Mean_statistics.xlsx','Sheet',model);
    writematrix(CI,'CI.xlsx','Sheet',model);
    writematrix(Daily_timestep_mean,'Daily_timestep_mean.xlsx','Sheet',model);
    writetable(w_listLL,'LL.xlsx','Sheet',model);
    writetable(AIC_final_list,'AIC.xlsx','Sheet',model);
    writetable(BIC_final_list,'BIC.xlsx','Sheet',model);
    writetable(KGE_final_list,'KGE.xlsx','Sheet',model);
    writetable(KGEi_final_list,'KGEi.xlsx','Sheet',model);
    writetable(KGEm_final_list,'KGEm.xlsx','Sheet',model);
    writetable(NSE_final_list,'NSE.xlsx','Sheet',model);
    writetable(NSEi_final_list,'NSEi.xlsx','Sheet',model);
    writetable(NSElog_final_list,'NSElog.xlsx','Sheet',model);
    writetable(RMSE_final_list,'RMSE.xlsx','Sheet',model);
    writetable(NRMSE_final_list,'NRMSE.xlsx','Sheet',model);
    writetable(R_final_list,'R.xlsx','Sheet',model);

    clear w_listLL AIC_final_list BIC_final_list KGE_final_list KGEi_final_list
    clear KGEm_final_list NSE_final_list NSEi_final_list NSElog_final_list
    clear RMSE_final_list NRMSE_final_list R_final_list Daily_timestep_mean
    clear b CI LL aic bic KGE KGEi KGEm NSEf NSEi NSElog RMSEE NRMSEE R qsim
    clear Par Q Q_sorted simLL simaic simbic simKGE simKGEi simKGEm simNSEf
    clear simNSEi simNSElog simRMSEE simNRMSEE simR
end

% MC graph
for i = 1:length(model_list)
    model = model_list{i};
    % Compare simulated and observed streamflow
    figure('color','w');
    box on;
    hold on;

    Q = Q_sim{i}.';
    for j = 1:numSample
        h(1) = plot(t,Q(:,j),'color',[0.5,0.5,0.5]);
    end
    h(2) = plot(t,input_climate.q,'r:');

    legend(h,'Simulated','Observed')
    title('Monte Carlo Simulations for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(h(2),'LineWidth',2)
    set(gca,'fontsize',16);
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename_MC = model_list{i}+"_MC.fig";
    savefig(gcf,filename_MC);
    close
    clear h ax Q
end

% Daily 90% uncertainty bounds of mean values
for i = 1:length(model_list)
    model = model_list{i};

    figure('color','w');
    box on;
    hold on;

    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    for j = 1:length(t)
        CI(j,1) = (Q_sorted(fifth_percentile,j));
        CI(j,2) = (Q_sorted(nintyfifth_percentile,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end

    sim_mean = Daily_timestep_mean(:,1);                                  % Daily Time-step Median of each model
    m(1) = area(t,CI(:,2),'FaceColor','blue','LineStyle','none','FaceAlpha',0.5);
    m(2) = area(t,CI(:,1),'FaceColor','white','LineStyle','none');

    % Plot daily time-step median
    m(3) = plot(t,sim_mean(:,1),'blue','LineWidth',1);
    m(4) = plot(t,input_climate.q,'r','LineWidth',1);
    %     sz = 10;
    %     m(4) = scatter(t,input_climate.q,sz,'red','filled');

    %m(1).HandleVisibility = 'off';
    m(2).HandleVisibility = 'off';

    legend('90% Uncertainty Bounds','Simulated','Observed')
    title('Q_realisations Mean for ',model,'Interpreter','none','FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
        'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
        '2009','2010','2011','2012','2013','2014','2015','2016'});
    set(gca,'fontsize',12);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_Q_realisations_mean = model_list{i}+"_Q_realisations_mean.fig";
    savefig(gcf,filename_Q_realisations_mean);                    % Export the Q_realisations graph as a figure

    close
    clear ax m Q Q_sorted CI
end

% Flow Duration Curve for mean series
for i = 1:length(model_list)
    model = model_list{i};
    Q_sorted = Q_sim{i};
    for j = 1:length(t)
        Q_sorted(:,j) = sort(Q_sorted(:,j));
    end

    Q = Q_sorted;
    clear Daily_timestep_mean
    for j = 1:length(t)
        Daily_timestep_mean(j,1) = mean(Q(fifth_percentile:nintyfifth_percentile,j));                           % Daily Time-step Mean of all the models
    end
    data_sim = Daily_timestep_mean.'; % change data here
    data_obs = (input_climate.q).';
    n_sim = length(data_sim); % count
    n_obs = length(data_obs); % count

    % Rank (order) for each value: same as rank function in EXCEL
    X = data_sim';
    Y = data_obs';
    for j = 1:length(X(1,:))
        [~, order_sim]= sort(X(:,j));
        m_sim(order_sim,j) = 1:length(X(:,j));
    end
    for j = 1:length(Y(1,:))
        [~, order_obs]= sort(Y(:,j));
        m_obs(order_obs,j) = 1:length(Y(:,j));
    end

    % rank is 'm_sim' and 'm_obs'
    % Exceedance Probability:
    for j=1:length(X)
        P_sim(j)=m_sim(j)./(n_sim+1);
    end
    for j=1:length(Y)
        P_obs(j)=m_obs(j)./(n_obs+1);
    end

    figure;
    box on;
    hold on;

    %Flow Duration Curve Plot
    l(1) = scatter(P_obs,Y,'red');
    l(2) = scatter(P_sim,X,'blue');

    legend(l,'Simulated','Observed')
    title('Flow Duration Curve for Ensemble Mean of ',model,'Interpreter','none','FontWeight','bold');
    xlabel('Exceedance Probability');
    ylabel('Discharge [mm]');
    set(gca,'fontsize',16);
    grid on

    ax = gca;
    ax.Layer = 'top';

    filename_FDC = model+"_FDC_mean.fig";
    savefig(gcf,filename_FDC);

    close
    clear l ax Q Q_sorted data_sim data_obs X Y m_obs m_sim P_sim P_obs
    clear n_sim n_obs filename_MC filename_FDC
end
