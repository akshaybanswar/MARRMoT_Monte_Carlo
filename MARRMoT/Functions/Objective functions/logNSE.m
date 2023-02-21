% logNSE Funktion

function lognse = logNSE( obs,sim )

%% invert the time series and add a small constant to avoid issues with 0 flows
% Pushpalatha et al (2012) suggests to set e at 1/100th of the mean of the
% observed flow, which is what we'll follow here. The constant is added
% before transforming flows.

% Find the constant
e = mean(obs)/100;

% Apply the constant and transform flows
obs1 = (obs+e);
sim1 = (sim+e);
%%
up=sum((log(sim1)-log(obs1)).^2);
down=sum((log(obs1)-log(mean(obs1))).^2);
lognse=1-(up/down);
end
