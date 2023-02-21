function [best_e,best_par,e_list] = efficiency_rank(efficiency,numSample,efficiency_name,parameters_set)

% INPUT
% 'efficiency' should be of 1 x n array
% 'numSample' should either be a number or a cell of 1 x n array
% 'parameters_set' should be a matrix of n x y array where y is the number
% of parameters and n is the number of samples
% 'efficiency_name' : 'AIC', 'BIC', 'LL' or 'LOGLIKELIHOOD', 'KGE', 'NSE',
%                     'R' or 'PEARSON', 'RMSE', 'NRMSE', 'GLUE'
%
% OUTPUT
% 'best_e' is the best value chosen among 1 x n array of efficiency
% 'best_par' returns the best parameter set corresponding to the 'best_e'
% 'e_list' gives the table from best to worst values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIC   : minimum of weighted AIC                                         %
% BIC   : minimum of BIC                                                  %
% LL    : minimum of LL                                                   %
% KGE   : -1 < KGE < 1                                                    %
% NSE   : -Inf < NSE < 1                                                  %
% R     : maximum of R                                                    %
% RMSE  : minimum of RMSE                                                 %
% NRMSE : minimum of NRMSE                                                %
% GLUE  : minimum of GLUE                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch efficiency_name

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'AIC'

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        min_aic = nanmin(efficiency);

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idx(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Calculate and check if deltaAIC value is less than 0
        for j = 1:length(idx)
            x = idx(j);
            delta_aic(1,j) = efficiency(1,x) - min_aic;
            par(j,:) = parameters_set(x,:);
            if delta_aic(1,j) < 0
                min_aic = efficiency(1,x);
                d = 1;
                for k = 1:j
                    y = idx(k);
                    delta_aic(1,k) = efficiency(1,y) - min_aic;
                    par(k,:) = parameters_set(y,:);
                    d = d + 1;
                end
            end
        end

        % Calculating the relative likelihood
        for j = 1:length(idx)
            idy(j) = j;
            relative_likelihood(1,j) = exp(-0.5*delta_aic(1,j));
        end

        relative_likelihood_sum = sum(relative_likelihood(1,:));

        for j = 1:length(idy)
            x = idy(j);
            w(j) = relative_likelihood(1,j)/relative_likelihood_sum;              % Calculating Akaike Weights
            weighted_e(j) = w(j) * efficiency(1,x);
        end

        max_w = max(weighted_e);

        for j = 1:length(idy)
            y = idy(j);
            if weighted_e(j) == max_w
                best_e = efficiency(1,y);
                best_par = par(j,:);
            end
        end

        e_sorted = weighted_e;
        w_sorted = w;
        par_sorted = par;

        for j = 1:length(idy)
            y = idy(j);
            ic_sorted(j) = efficiency(1,y);
        end

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i-1) < e_sorted(i))
                    e_sorted = swap(e_sorted,i-1,i);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    w_sorted = swap(w_sorted,i-1,i);
                    ic_sorted = swap(ic_sorted,i-1,i);
                    models_sorted = swap(models_sorted,i-1,i);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(w_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = ic_sorted(1,j);
            e_list{j,5} = w_sorted(1,j);
            e_list{j,6} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'Weighted AIC';
        e_list.Properties.VariableNames{4} = 'AIC';
        e_list.Properties.VariableNames{5} = 'Weights';
        e_list.Properties.VariableNames{6} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%%%% BIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'BIC'

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idy(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(1,j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
        end

        min_e = nanmin(e);
        for j = 1:length(idy)
            if e(j) == min_e
                y = idy(1,j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i) < e_sorted(i-1))
                    e_sorted = swap(e_sorted,i,i-1);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted = swap(models_sorted,i,i-1);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'BIC';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%% LL, Loglikelihood %%%%%%%%%%%%%%%%%%%%%%%%%%

    case {'LL' , 'LOGLIKELIHOOD'}

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idy(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
        end

        max_e = max(e);
        for j = 1:length(idy)
            if e(j) == max_e
                y = idy(j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i-1) < e_sorted(i))
                    e_sorted = swap(e_sorted,i-1,i);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted = swap(models_sorted,i-1,i);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'Loglikelihood';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% KGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'KGE'

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idx(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if KGE value lies between -0.41 to 1 and storing their location
        a = 1;
        b = 0;
        for j = 1:length(idx)
            x = idx(j);
            if efficiency(1,x) >= -0.41 && efficiency(1,x) <= 1
                idy(a) = idx(j);
                models_sorted1{a} = models_sorted{x};
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are not in the range
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
            delta(1,j) = 1 - e(j);
        end

        min_delta = min(delta);
        for j = 1:length(idy)
            if delta(1,j) == min_delta
                y = idy(j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i-1) < e_sorted(i))
                    e_sorted = swap(e_sorted,i-1,i);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted1 = swap(models_sorted1,i-1,i);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted1{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'KGE';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% NSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'NSE'

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idx(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if NSE value lies between -Inf to 1 and storing their location
        a = 1;
        c = 0;
        for j = 1:length(idx)
            x = idx(j);
            if efficiency(1,x) <= 1
                idy(a) = idx(j);
                models_sorted1{a} = models_sorted{x};
                a = a + 1;
            else
                c = c + 1;
            end
        end

        % Return NaN if all the values are not in the range
        if c == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
            delta(1,j) = 1 - e(j);
        end

        min_delta = min(delta);
        for j = 1:length(idy)
            if delta(1,j) == min_delta
                y = idy(j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i-1) < e_sorted(i))
                    e_sorted = swap(e_sorted,i-1,i);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted1 = swap(models_sorted1,i-1,i);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted1{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'NSE';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% RMSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'RMSE'

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idy(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(1,j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
        end

        min_e = nanmin(e);
        for j = 1:length(idy)
            if e(j) == min_e
                y = idy(1,j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i) < e_sorted(i-1))
                    e_sorted = swap(e_sorted,i,i-1);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted = swap(models_sorted,i,i-1);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'RMSE';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%%%%% NRMSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'NRMSE'

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idy(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(1,j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
        end

        min_e = nanmin(e);
        for j = 1:length(idy)
            if e(j) == min_e
                y = idy(1,j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i) < e_sorted(i-1))
                    e_sorted = swap(e_sorted,i,i-1);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted = swap(models_sorted,i,i-1);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'NRMSE';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%% R, PEARSON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {'R' , 'PEARSON'}

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idx(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if R value lies between -1 to 1 and storing their location
        a = 1;
        c = 0;
        for j = 1:length(idx)
            x = idx(j);
            if efficiency(1,x) >= -1 && efficiency(1,x) <= 1
                idy(a) = idx(j);
                models_sorted1{a} = models_sorted{x};
                a = a + 1;
            else
                c = c + 1;
            end
        end

        % Return NaN if all the values are not in the range
        if c == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
        end

        max_e = max(e);
        for j = 1:length(idy)
            if e(j) == max_e
                y = idy(j);
                best_e = e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i-1) < e_sorted(i))
                    e_sorted = swap(e_sorted,i-1,i);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    models_sorted1 = swap(models_sorted1,i-1,i);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted1{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'R';
        e_list.Properties.VariableNames{4} = 'Parameter Set';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {'GLUE'}

        if ~exist('parameters_set','var')
            % third parameter does not exist, so default it to something
            parameters_set(1:numSample,1) = NaN;
        end

        if ~isreal(efficiency)
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        % Check if the numSample passed is a number or a cell array of the list of models
        y = isnumeric(numSample);
        if y == 1
            nSample = numSample;
            for j = 1:nSample
                models_sorted{j} = j;
            end
        else
            z = iscell(numSample);
            if z == 1
                nSample = length(numSample);
                models_sorted = numSample;
            else
                disp('Number of Samples should either be a number or a cell array.')
                return
            end
        end

        % Filter out the NaN values
        a = 1;
        b = 0;
        for j = 1:nSample
            q = isnan(efficiency(1,j));
            if q == 0
                idy(a) = j;
                a = a + 1;
            else
                b = b + 1;
            end
        end

        % Return NaN if all the values are NaN
        if b == j
            best_e = NaN;
            best_par = NaN;
            e_list = NaN;
            e_list = array2table(e_list);
            return
        end

        for j = 1:length(idy)
            y = idy(j);
            e(j) = efficiency(1,y);
            par(j,:) = parameters_set(y,:);
        end

        sum_e = sum(e);
        max_e = max(e);

        for j = 1:length(idy)
            x = idy(j);
            w(j) = e(j)/sum_e;
            weighted_e(j) = w(j)*e(x);
        end

        max_weighted_e = max(weighted_e);

        for j = 1:length(idy)
            if weighted_e(j) == max_weighted_e
                y = idy(j);
                best_e = weighted_e(j);
                best_par = par(j,:);
            end
        end

        e_sorted = weighted_e;
        w_sorted = w;
        LL = e;
        par_sorted = par;

        % Bubble sort
        n = length(e_sorted);
        while (n > 0)
            % Iterate through x
            nnew = 0;
            for i = 2:n
                % Swap elements in wrong order
                if (e_sorted(i-1) < e_sorted(i))
                    e_sorted = swap(e_sorted,i-1,i);
                    val = par_sorted(i,:);
                    par_sorted(i,:) = par_sorted(i-1,:);
                    par_sorted(i-1,:) = val;
                    w_sorted = swap(w_sorted,i-1,i);
                    LL = swap(LL,i-1,i);
                    models_sorted = swap(models_sorted,i-1,i);
                    nnew = i;
                end
            end
            n = nnew;
        end

        for j = 1:length(e_sorted)
            e_list{j,1} = j;
            e_list{j,2} = models_sorted{j};
            e_list{j,3} = e_sorted(j);
            e_list{j,4} = LL(j);
            e_list{j,5} = w_sorted(j);
            e_list{j,6} = par_sorted(j,:);
        end

        e_list = cell2table(e_list);
        e_list.Properties.VariableNames{1} = 'Rank';
        e_list.Properties.VariableNames{2} = 'Model Name';
        e_list.Properties.VariableNames{3} = 'Weighted Loglikelihood';
        e_list.Properties.VariableNames{4} = 'Loglikelihood';
        e_list.Properties.VariableNames{5} = 'GLUE Weights';
        e_list.Properties.VariableNames{6} = 'Parameter Set';

end
end

function x = swap(x,i,j)
% Swap x(i) and x(j)
% Note: In practice, x should be passed by reference
val = x(i);
x(i) = x(j);
x(j) = val;
end

