% This file is part of the Modular Assessment of Rainfall-Runoff Models 
% Toolbox (MARRMoT) – User manual. It contains an example application of a
% single model to a single catchment. See section 3 in the User Manual 
% for details.
% 
% Author:   Wouter J.M. Knoben
% Date:     26-09-2018
% Contact:  w.j.m.knoben@bristol.ac.uk
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% This example workflow includes 5 steps:
%
% 1. Data preparation
% 2. Model choice and setup
% 3. Model solver settings and time-stepping scheme
% 4. Model runs
% 5. Output vizualization

%% 1. Prepare data
% Load the data
load fluxInput.mat
load Q_Tua.mat

% Create a climatology data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_climatology.precip   = fluxInput.precip;      % Daily data: P rate  [mm/d]
input_climatology.temp     = fluxInput.temp;        % Daily data: mean T  [degree C]
input_climatology.pet      = fluxInput.pet;         % Daily data: Ep rate [mm/d]
input_climatology.q        = Q_Tua.discharge;       % Daily data: Run off [mm/d]
input_climatology.delta_t  = 1;                     % time step size of the inputs: 1 [d]

%% 2. Define the model settings
% NOTE: this example assumes that parameter values for this combination of
% model and catchment are known. 

% Model name 
% NOTE: these can be found in the Model Descriptions
model = 'm_07_gr4j_4p_2s';                     

% Parameter values
% NOTE: descriptions of these parameters can be found in the Model
% descriptions (supplementary materials to the main paper). Alternatively,
% the parameters are described in each model function. Right-click the
% model name above (i.e. 'm_29_hymod_5p_5s') and click 
% "Open 'm_29_hymod_5p_5s'". Parameters are listed on lines 44-50.

input_theta     = [ 1.0;                                                    
                     -7.4;                                                   
                     136.4;                                                                                        
                     2.0];                                                 

% Initial storage values
% NOTE: see the model function for the order in which stores are given. For
% HyMOD, this is on lines 86-91.

input_s0       = [15;                                                       % Initial soil moisture storage [mm]
                  10];                                                      % Initial slow flow storage [mm]

%% 3. Define the solver settings  
% Create a solver settings data input structure. 
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver.name              = 'createOdeApprox_IE';                      % Use Implicit Euler to approximate ODE's
input_solver.resnorm_tolerance = 0.1;                                       % Root-finding convergence tolerance
input_solver.resnorm_maxiter   = 6;                                         % Maximum number of re-runs


%% 4. Run the model and extract all outputs
% This process takes ~6 seconds on a i7-4790 CPU 3.60GHz, 4 core processor.
% NOTE: requesting a water balance check by using the fourth output
% argument automatically prints a water balance overview to the screen.
[output_ex,...                                                              % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
 output_in,...                                                              % Internal model fluxes
 output_ss,....                                                             % Internal storages
 output_waterbalance] = ...                                                 % Water balance check
                    feval(model,...                                         % Model function name
                          input_climatology,...                             % Time series of climatic fluxes in simulation period
                          input_s0,...                                      % Initial storages
                          input_theta,...                                   % Parameter values
                          input_solver);                                    % Details of numerical time-stepping scheme
    
%% 5. Analyze the outputs                   
% Prepare a time vector
t = datetime(1972,01,01):datetime(2016,12,31);

% Compare simulated and observed streamflow by calculating the Kling-Gupta
% Efficiency (KGE). Other objective functions provided are inverse KGE and
% multi-objective average KGE (0.5*(KGE(Q) + KGE(1/Q))
tmp_obs  = input_climatology.q;
tmp_sim  = output_ex.Q;
tmp_kge  = of_KGE(tmp_obs,tmp_sim);                                         % KGE on regular flows
tmp_kgei = of_inverse_KGE(tmp_obs,tmp_sim);                                 % KGE on inverse flows
tmp_kgem = of_mean_hilo_KGE(tmp_obs,tmp_sim);                               % Average of KGE(Q) and KGE(1/Q)

KGEi_cal = of_inverse_KGE(tmp_obs(367:9132,1),tmp_sim(1,367:9132).');
KGEi_eval = of_inverse_KGE(tmp_obs(9499:16437,1),tmp_sim(1,9499:16437).');

figure('color','w'); 
    box on;
    hold on; 
    
    h1 = plot(t,tmp_obs);
    h2 = plot(t,tmp_sim);
    
    legend('Observed','Simulated')
    title(['Kling-Gupta Efficiency = ',num2str(tmp_kge)])
    ylabel('Streamflow [mm/d]')
    xlabel('Time [d]')
    datetick;
    set(h1,'LineWidth',2)
    set(h2,'LineWidth',2)
    set(gca,'fontsize',16);

clear hi h2 
    
% Investigate internal storage changes
figure('color','w');
    
    p1 = subplot(311);
        hold on;
        h1 = plot(t,output_ss.S1);
        title('Simulated storages')
        ylabel('Soil moisture [mm]')
        datetick;
        
    p2 = subplot(312);
        box on;
        hold all;
        h2 = plot(t,output_ss.S2);
        h3 = plot(t,output_ss.S3);
        h4 = plot(t,output_ss.S4);
        legend('Fast store 1','Fast store 2','Fast store 3')
        ylabel('Fast stores [mm]')
        datetick;
        
    p3 = subplot(313);
        h5 = plot(t,output_ss.S5);
        ylabel('Slow store [mm]')
        xlabel('Time [d]')
        datetick;

    set(p1,'fontsize',16)
    set(p2,'fontsize',16)
    set(p3,'fontsize',16)
        
    set(h1,'LineWidth',2)
    set(h2,'LineWidth',2)
    set(h3,'LineWidth',2)
    set(h4,'LineWidth',2)
    set(h5,'LineWidth',2)

clear p* h* t










