%% The workflow includes 5 steps:
%
% 1. Data preparation
% 2. Model choice and setup
% 3. Model solver settings and time-stepping scheme
% 4. Monte Carlo run
% 5. Extracting and formating the output in the desired format
% 6. Output vizualization

tic

%% 1. Prepare data
clear all; close all; clc;
%addpath(genpath('C:\Eddys\MonteCarlo\MARRMoT_Monte_Carlo'));

% Load the data
load INPUT.mat

% Create a climatology data input structure.
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_climatology.precip   = INPUT.p;                                      % Daily data: P rate  [mm/d]
input_climatology.temp     = INPUT.temp;                                   % Daily data: mean T  [degree C]
input_climatology.pet      = INPUT.pet;                                    % Daily data: Ep rate [mm/d]
input_climatology.q        = (INPUT.qobs*86400)/(INPUT.area*1000);         % Daily data: Run off [mm/d]
input_climatology.delta_t  = 1;                                            % Time step size of the inputs: 1 [d]
input_climatology.date     = INPUT.date;                                   % Prepare a time vector

% Time periods for calibration and evaluation.
% Note: generally a 'warm-up period' is used to lessen the impact of the
% initial conditions. Examples include running 1 year of data iteratively
% until model stores reach an equilibrium, or choosing an arbitrary cut-off
% point before which the simulations are judged to be inaccurate and only
% after which the objective function is calculated. For brevity, this step
% is ignored in this example and the full calibration and evaluation
% periods are used to calculate the objective function. Initial storages
% are estimated as 0 (see line 60 of this script).
time_start      = 1;                                                       % 01.01.2000 (730486)
time_end        = 6210;                                                    % 31.12.2016 (736695)
warmup          = 366;                                                     % NUMBER OF TIME STEPS TO DISREGARD IN THE COMPUTATION OF THE FITNESS CRITERTA #WOEHLING#

% Create temporary calibration time series
input_climate.precip    = input_climatology.precip(time_start:time_end);
input_climate.temp      = input_climatology.temp(time_start:time_end);
input_climate.pet       = input_climatology.pet(time_start:time_end);
input_climate.delta_t   = input_climatology.delta_t;
input_climate.q         = input_climatology.q(time_start:time_end);
input_climate.date      = input_climatology.date;
qobs                    = input_climate.q;
numObs                  = length(find(input_climatology.q(warmup:end) >= 0));

%% 2. Define the model settings and number of samples
% NOTE: this example assumes that the model parameters for each model will
% be sampled as part of the investigation. See lines 69-77 in this script.

% Model name
% NOTE: these can be found in the Model Descriptions.
model_list  = {'m_05_ihacres_7p_1s',...
    'm_07_gr4j_4p_2s',...
    'm_14_topmodel_7p_2s',...
    'm_26_flexi_10p_4s',...
    'm_35_mopex5_12p_5s',...
    'm_36_modhydrolog_15p_5s',...
    'm_37_hbv_15p_5s',...
    'm_44_echo_16p_6s'};

% Define the requested number of samples
numSample = 10000;

%% 3. Define the solver settings
% Create a solver settings data input structure.
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver.name              = 'createOdeApprox_IE';                     % Use Implicit Euler to approximate ODE's
input_solver.resnorm_tolerance = 0.1;                                      % Root-finding convergence tolerance
input_solver.resnorm_maxiter   = 6;                                        % Maximum number of re-runs

%% 4. Get Monte Carlo Sampling Results
% Run Monte Carlo Analysis
% NOTE: Using Latin Hypercube Sampling
% Start the sampling
% Prepare a time vector
t = INPUT.date;

% Preallocation of variables
for i = 1:length(model_list)
    %input_parameter{i} = NaN;
    for j = 1:numSample
        QEA{1,i}{1,j}.Q(1:length(t)) = NaN;
        QEA{1,i}{1,j}.Ea(1:length(t)) = NaN;
    end
end

for i = 1:length(model_list)

    % Display progress
    disp(['Now starting model ',model_list{i},'.'])
    disp(' ')

    % Define the current model
    model = model_list{i};

    [input_parameter{i},...                                                % Parameters set using Latin Hypercube Sampling
        QEA{i},...                                                         % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
        ~,...                                                              % Internal model fluxes
        ~,...                                                              % Internal storages
        ~] = ...                                                           % Water balance check
        MonteCarloSamplingResults(numSample, ...                           % Define the requested number of samples
        model, ...                                                         % List of model names
        input_climate, ...                                                 % Time series of climatic fluxes in simulation period
        input_solver);                                                     % Solver Settings
end

% SAVE RESULTS
dumpall = ['save everything_' num2str(now) '.mat'];
eval(dumpall);

%% 5. Extracting and Formating the Desired Output
% Extracting the discharge and parameter values from the results
% Formating the values to fit in a tabular format to generate excel files for further use

for i = 1:length(model_list)

    Q(1:numSample,1:size(input_climatology.precip,1)) = NaN;
    for j=1:numSample
        Q(j,:)=QEA{1,i}{1,j}.Q;
    end
    %Par = cell2mat(input_parameter{1,i});
    Q_sim{i} = Q;
end
fifth_percentile = round(((5/100)*numSample))+1;
nintyfifth_percentile = round(((95/100)*numSample))-1;

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
    model_range = feval([model,'_parameter_ranges']);                          % Call the function that stores parameter ranges

    % Find number of stores and parameters
    numPar    = size(model_range,1);

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
        RMSEE(1,j)    = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
        NRMSEE(1,j)   = RMSEE(1,j)/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
        R(1,j)        = corr(qobs(warmup:end,1),qsim(1,warmup:end).');
    end
    Par = (input_parameter{1,i});
    [~,~,w_listLL]             = efficiency_rank(LL(1,:),numSample,'LL',Par);
    [~,~,AIC_final_list]       = efficiency_rank(aic(1,:),numSample,'AIC',Par);
    [~,~,BIC_final_list]       = efficiency_rank(bic(1,:),numSample,'BIC',Par);
    [~,~,KGE_final_list]       = efficiency_rank(KGE(1,:),numSample,'KGE',Par);
    [~,~,KGEi_final_list]      = efficiency_rank(KGEi(1,:),numSample,'KGE',Par);
    [~,~,KGEm_final_list]      = efficiency_rank(KGEm(1,:),numSample,'KGE',Par);
    [~,~,NSE_final_list]       = efficiency_rank(NSEf(1,:),numSample,'NSE',Par);
    [~,~,NSEi_final_list]      = efficiency_rank(NSEi(1,:),numSample,'NSE',Par);
    [~,~,NSElog_final_list]    = efficiency_rank(NSElog(1,:),numSample,'NSE',Par);
    [~,~,RMSE_final_list]      = efficiency_rank(RMSEE(1,:),numSample,'RMSE',Par);
    [~,~,NRMSE_final_list]     = efficiency_rank(NRMSEE(1,:),numSample,'NRMSE',Par);
    [~,~,R_final_list]         = efficiency_rank(R(1,:),numSample,'R',Par);

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
    simRMSEE    = rmse(qobs(warmup:end,1),qsim(1,warmup:end).');
    simNRMSEE   = simRMSEE/(max(qobs(warmup:end,1)) - min(qobs(warmup:end,1)));
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

%% 6. Analyse the Output
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

toc