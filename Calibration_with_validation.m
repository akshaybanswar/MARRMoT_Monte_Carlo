% This file is part of the Modular Assessment of Rainfall-Runoff Models
% Toolbox (MARRMoT) – User manual. It contains an example application of
% calibration of a single models to a single catchment. See section 3 in
% the User Manual for details.
%
% NOTE: this example uses a custom function 'fminsearchbnd', which is a
% basic optimization that lets the user specify constraints in the solution
% space. In this example these constraints are the model's parameter
% ranges. The file can be downloaded from Matlab's File Exchange:
% https://uk.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
%
% NOTE: this file does not work very well in Octave. Octave users might
% need to consider alternative parameter optimisation methods.
%
% Author:   Wouter J.M. Knoben
% Date:     26-09-2018
% Contact:  w.j.m.knoben@bristol.ac.uk
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% This example workflow includes 8 steps:
%
% 1. Check whether the optimization function is present
% 2. Data preparation
% 3. Model choice and setup
% 4. Model solver settings and time-stepping scheme
% 5. Calibration settings
% 6. Model calibration
% 7. Evaluation of calibration results
% 8. Output vizualization
clear
tic
%% 1. Check function requirements
% Check whether fminsearchbnd is properly installed
if ~exist('fminsearchbnd','file')
    disp(['Optimizer ''fminsearchbnd'' not found. Visit the following ',...
        'link to download: ',...
        'https://uk.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon']);
    return
end

%% 2. Prepare data
% Load the data
load INPUT.mat

% Create a climatology data input structure.
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_climatology.precip   = INPUT.p;                                       % Daily data: P rate  [mm/d]
input_climatology.temp     = INPUT.temp;                                    % Daily data: mean T  [degree C]
input_climatology.pet      = INPUT.pet;                                     % Daily data: Ep rate [mm/d]
input_climatology.q        = (INPUT.qobs*86400)/(INPUT.area*1000);          % Daily data: Run off [mm/d]
input_climatology.delta_t  = 1;                                             % time step size of the inputs: 1 [d]

%% 3. Define the model settings
model_list  = {'m_05_ihacres_7p_1s',...
    'm_07_gr4j_4p_2s',...
    'm_14_topmodel_7p_2s',...
    'm_26_flexi_10p_4s',...
    'm_35_mopex5_12p_5s',...
    'm_36_modhydrolog_15p_5s',...
    'm_37_hbv_15p_5s',...
    'm_44_echo_16p_6s'};

%% 4. Define the solver settings
% Create a solver settings data input structure.
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver.name              = 'createOdeApprox_IE';                      % Use Implicit Euler to approximate ODE's
input_solver.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
input_solver.resnorm_maxiter   = 6;                                         % Maximum number of re-runs

%% 5. Define calibration settings
% Settings for 'fminsearchbnd'
optim_settings = optimset(...                                               % Using 'optimset' ensures that all non-specified fields are included in the resulting settings structure
    'Display','iter',...                                    % Track progress (default = 'off')
    'MaxIter',1000,...                                        % Stop after this many iterations (default = 200*numPar)
    'MaxFunEvals',1000,...                                   % Stop after this many function evaluations (default = 200*numPar)
    'TolFun',1e-2,...                                       % Stop if objective function change is below this tolerance (default = 1e-4)
    'TolX',1e-2);                                           % Stop if changes in parameter values is below this tolerance (default = 1e-4)

% Choose the objective function
of_name      = 'of_KGE';                                                    % This function is provided as part of MARRMoT. See ./MARRMoT/Functions/Objective functions

% Time periods for calibration and evaluation.
% Note: generally a 'warm-up period' is used to lessen the impact of the
% initial conditions. Examples include running 1 year of data iteratively
% until model stores reach an equilibrium, or choosing an arbitrary cut-off
% point before which the simulations are judged to be inaccurate and only
% after which the objective function is calculated. For brevity, this step
% is ignored in this example and the full calibration and evaluation
% periods are used to calculate the objective function. Initial storages
% are estimated as 0 (see line 60 of this script).
time_cal_start  = 1;                   % 01.01.2000 (730486)
%time_cal_end    = 6210;                % 31.12.2016 (736695)
time_cal_end    = 3105;                % 01.06.2008 (733590)
time_eval_start = 3106;                % 02.06.2008 (733591)
time_eval_end   = 6210;                % 31.12.2016 (736695)
warmup          = 366;              % NUMBER OF TIME STEPS TO DISREGARD IN THE COMPUTATION OF THE FITNESS CRITERTA #WOEHLING#

%% 6. Calibrate the model
% Each MARRMoT model provides its outputs in a standardized way. We need
% the simulated flows for a given parameter set, so that we can compare
% simulated and observed flow, and determine the fitness of a particular
% parameter set. Simulated flow is stored in a structure of the form
% output.Q. The auxiliary function 'workflow_calibrationAssist' takes our
% specified model, inputs and a parameter set, runs the model and returns
% just simulated flow. The returned variable is used to optimize the
% parameter values. See lines 113-142.

% Create temporary calibration time series
input_climate_cal.precip  = input_climatology.precip(time_cal_start:time_cal_end);
input_climate_cal.temp    = input_climatology.temp(time_cal_start:time_cal_end);
input_climate_cal.pet     = input_climatology.pet(time_cal_start:time_cal_end);
input_climate_cal.delta_t = input_climatology.delta_t;
q_obs_cal                 = input_climatology.q(time_cal_start:time_cal_end);

for i = 1:length(model_list)
    model    = model_list{i};                                               % Name of the model function (these can be found in Supporting Material 2)
    parRange = feval([model,'_parameter_ranges']);                              % Parameter ranges
    numPar   = size(parRange,1);                                                % Number of parameters
    numStore = str2double(model(end-1));                                        % Number of stores
    input_s0 = zeros(numStore,1);                                               % Initial storages (see note in paragraph 5 on model warm-up)

    % Create a temporary function that calls the model, using the assisting
    % function so that we get Qsim as an output, having only the parameter
    % values as unknown input
    q_sim_fun = @(par) workflow_calibrationAssist(...                           % Auxiliary function that returns Qsim
        model,...                                            % Model name
        input_climate_cal,...                                % Climate data during calibration period
        input_s0,...                                         % Initial storages
        par,...                                              % We want to optimize these
        input_solver);                                       % Solver settings

    % Create the objective function we want to optimize. 'fminsearchbnd' is a
    % minimizer and the objective function 'of_KGE' has range <-Inf,1] (<bad
    % performance, good performance]. We thus want to optimize (minimize)
    % -1*of_KGE which has range [-1,Inf> ([good,bad>).
    cal_fun = @(par) -1*(...
        feval(of_name, ...                                         % The objective function, here this is 'of_KGE'
        q_obs_cal,...                                        % Observed flow during calibration period
        q_sim_fun(par), ...                                  % Simulated flow for a given parameter set
        warmup));

    % Create initial guesses for the optimizer(LatinHypercube Sampling)
    par_ini = parRange(:,1)+lhsdesign(numPar,1).*(parRange(:,2)-parRange(:,1));

    % Call 'fminsearchbound' to optimize parameters
    disp('--- Calibration starting ---')
    [par_opt,of_cal] = fminsearchbnd(...                                        % Returns optimized parameters and the objective function value
        cal_fun,...                                             % Function that compares Qobs and Qsim for a given parameter set
        par_ini,...                                             % Initial parameter guess
        parRange(:,1),...                                       % Lower parameter bounds
        parRange(:,2),...                                       % Upper parameter bounds
        optim_settings);                                        % Settings for fminsearchbnd

    %% 7. Evaluate the calibrated parameters on unseen data
    % Create temporary evaluation time series
    input_climate_eval.precip  = input_climatology.precip(time_eval_start:time_eval_end);
    input_climate_eval.temp    = input_climatology.temp(time_eval_start:time_eval_end);
    input_climate_eval.pet     = input_climatology.pet(time_eval_start:time_eval_end);
    input_climate_eval.delta_t = input_climatology.delta_t;
    q_obs_eval                 = input_climatology.q(time_eval_start:time_eval_end);

    % Run the model with calibration parameters
    model_out_eval = feval(model,...
        input_climate_eval,...                                    % Climate data during evaluation period
        input_s0,...                                              % Initial storages
        par_opt,...                                               % Calibrated parameters
        input_solver);                                            % Solver settings

    % Compute evaluation performance
    of_eval = feval(of_name,...                                                 % Objective function name (here 'of_KGE')
        q_obs_eval,...                                              % Observed flow during evaluation period
        model_out_eval.Q, ...                                       % Simulated flow during evaluation period, using calibrated parameters
        warmup);

    %% 8. Visualise the results
    % Get simulated flow during calibration
    model_out_cal = feval(model,...
        input_climate_cal,...                                     % Climate data during evaluation period
        input_s0,...                                              % Initial storages
        par_opt,...                                               % Calibrated parameters
        input_solver);                                            % Solver settings

    %% #WOEHLING#: COMPUTE SEVERAL FITNESS CRITERIA UNDER CONSIDERATION OF A SPIN-UP PHASE
    %     tmp_obs  = input_climatology.q;
    KGEi_cal  = of_inverse_KGE(q_obs_cal,model_out_cal.Q(1,:).',warmup);
    KGE_cal  = of_KGE(q_obs_cal,model_out_cal.Q(1,:).',warmup);
    KGEm_cal = of_mean_hilo_KGE(q_obs_cal,model_out_cal.Q(1,:).',warmup);
    NSE_cal  = NSE(q_obs_cal(warmup:end,1),model_out_cal.Q(1,warmup:end).');
    RMSE_cal  = rmse(q_obs_cal(warmup:end,1),model_out_cal.Q(1,warmup:end).');
    NRMSE_cal  = RMSE_cal/(max(q_obs_cal(warmup:end,1)) - min(q_obs_cal(warmup:end,1)));
    R_cal  = corr(q_obs_cal(warmup:end,1),model_out_cal.Q(1,warmup:end).');
    NSElog_cal  = logNSE(q_obs_cal(warmup:end,1),model_out_cal.Q(1,warmup:end).');
    NSEi_cal  = inverse_NSE(q_obs_cal(warmup:end,1),model_out_cal.Q(1,warmup:end).');
    AIC_cal = of_aic(q_obs_cal,model_out_cal.Q(1,:).',warmup,numPar);
    BIC_cal = of_bic(q_obs_cal,model_out_cal.Q(1,:).',warmup,numPar);
    LL_cal = of_loglikelihood(q_obs_cal,model_out_cal.Q(1,:).',warmup,numPar);

    KGEi_eval  = of_inverse_KGE(q_obs_eval,model_out_eval.Q(1,:).',warmup);
    KGE_eval  = of_KGE(q_obs_eval,model_out_eval.Q(1,:).',warmup);
    KGEm_eval = of_mean_hilo_KGE(q_obs_eval,model_out_eval.Q(1,:).',warmup);
    NSE_eval  = NSE(q_obs_eval(warmup:end,1),model_out_eval.Q(1,warmup:end).');
    RMSE_eval  = rmse(q_obs_eval(warmup:end,1),model_out_eval.Q(1,warmup:end).');
    NRMSE_eval  = RMSE_eval/(max(q_obs_eval(warmup:end,1)) - min(q_obs_eval(warmup:end,1)));
    R_eval  = corr(q_obs_eval(warmup:end,1),model_out_eval.Q(1,warmup:end).');
    NSElog_eval  = logNSE(q_obs_eval(warmup:end,1),model_out_eval.Q(1,warmup:end).');
    NSEi_eval  = inverse_NSE(q_obs_eval(warmup:end,1),model_out_eval.Q(1,warmup:end).');
    AIC_eval = of_aic(q_obs_eval,model_out_eval.Q(1,:).',warmup,numPar);
    BIC_eval = of_bic(q_obs_eval,model_out_eval.Q(1,:).',warmup,numPar);
    LL_eval = of_loglikelihood(q_obs_eval,model_out_eval.Q(1,:).',warmup,numPar);
    FIT = rows2vars(table(KGE_cal, KGE_eval, KGEi_cal,KGEi_eval, KGEm_cal, KGEm_eval, ...
        NSE_cal, NSE_eval, NSElog_cal, NSElog_eval,...
        NSEi_cal, NSEi_eval, RMSE_cal, RMSE_eval, NRMSE_cal, NRMSE_eval, R_cal, R_eval,...
        AIC_cal, AIC_eval, BIC_cal, BIC_eval, LL_cal, LL_eval));
    qsimcal = model_out_cal.Q(1,:).';
    qsimeval = model_out_eval.Q(1,:).';
    FIT.Properties.VariableNames{2} = model;

    writetable(FIT,"FIT.xlsx",'Sheet',model);
    writematrix(qsimcal,'qsim_cal.xlsx','Sheet',model);
    writematrix(qsimeval,'qsim_eval.xlsx','Sheet',model);

    %% Invert the calibration performance metric so that it once again has range
    % <-Inf,1]
    of_cal = -1*of_cal;

    % Prepare a time vector
    t = INPUT.date;

    % Compare simulated and observed
    figure('color','w');
    box on;
    hold all;

    % Flows
    l(1) = plot(t(time_cal_start:time_cal_end),q_obs_cal,'k');              % Qobs_cal
    plot(t(time_eval_start:time_eval_end),q_obs_eval,'k')                   % Qobs_eval
    l(2) = plot(t(time_cal_start:time_cal_end),model_out_cal.Q,'r');        % Qsim_cal
    plot(t(time_eval_start:time_eval_end),model_out_eval.Q,'r')             % Qsim_eval

    % Dividing line
    l(3) = plot([t(time_cal_end),t(time_cal_end)],[0,170],'--b','linewidth',2);

    % Legend & text
    l = legend(l,'Observed', 'Simulated','Cal // Eval');
    title('Model calibration results for ', model,'FontWeight','bold')
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')

    txt_cal  = sprintf('Calibration period \nKGE = %.2f ',of_cal);
    txt_eval = sprintf('Evaluation period \nKGE = %.2f',of_eval);
    text(t(time_cal_end-450),57,txt_cal,'fontsize',12);
    text(t(time_eval_start+310),57,txt_eval,'fontsize',12);
    set(gca,'fontsize',12);

    % Other settings
    datetick;
    xlim([730486 736695]);
    set(gca,'fontsize',12,'xtick',[730486 730852 731217 731582 731947 732313 732678 ...
        733043 733408 733774 734139 734473 734868 735235 735600 735965 736330],...
       'xticklabel',{'2000','2001','2002','2003','2004','2005','2006','2007','2008',...
       '2009','2010','2011','2012','2013','2014','2015','2016'});
    ylim([0,60])
    set(gca,'TickLength',[0.005,0.005])
    grid on;
    ax = gca;
    ax.Layer = 'top';

    filename2 = model_list{i}+"_Calibration.fig";
    savefig(gcf,filename2);

    clear ax l


end

toc



