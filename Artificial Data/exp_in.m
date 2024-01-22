function [f, LL] = exp_in(lam, data, window)
% Likelihood contribution of the exponential distribution for incomplete 
% (letf-truncated and right-censored data), 

% data: the data to be fitting with 5 columns, ID, y1,y2,y3 and tag
% lam: lambda of exponential (log) 
% window: the window to be fitted, can be W1, W2 or W

%================================ preprocess initial parameters ============  
exponential_mean = 1 / exp(lam);
LL = [];   % loglikelihood of each gamma
for tag=1:8
    ind = data(:,5) == tag;  % indices to elements in 5th column of data that satisfy the equality
    temp_data = data(ind,:);
    [temp_n, ~] = size(temp_data); 
    
    % calculate pdf and survival function according to the 8 scenarios
    L = [];  % initial L
    if temp_n ~= 0
        y1 = temp_data(:,2);
        y2 = temp_data(:,3);
        y3 = temp_data(:,4);

        if window == "W1"
            if tag == 1
                term1 = exppdf(y1 + y2, exponential_mean);
                term2 = 1 - expcdf(y1, exponential_mean);
            elseif tag == 2
                term1 = exppdf(y2, exponential_mean);
            elseif tag == 3
                term1 = 1 - expcdf(y1 + y2, exponential_mean);
                term2 = 1 - expcdf(y1, exponential_mean);
            elseif tag == 4
                term1 = 1 - expcdf(y2, exponential_mean);
            elseif tag == 5
                term1 = 1 - expcdf(y2, exponential_mean);
            elseif tag == 6
                term1 = 1 - expcdf(y1 + y2, exponential_mean);
                term2 = 1 - expcdf(y1, exponential_mean);
            elseif tag == 7
                term1 = ones(temp_n,1);
            elseif tag == 8
                term1 = ones(temp_n,1);
            end
        end
        if window == "W2"
            if tag == 1
                term1 = ones(temp_n,1);
            elseif tag == 2
                term1 = ones(temp_n,1);
            elseif tag == 3
                term1 = exppdf(y1 + y2 + y3, exponential_mean);
                term2 = 1 - expcdf(y1 + y2, exponential_mean);
            elseif tag == 4
                term1 = exppdf(y2 + y3, exponential_mean);
                term2 = 1 - expcdf(y2, exponential_mean);
            elseif tag == 5
                term1 = 1 - expcdf(y2 + y3, exponential_mean);
                term2 = 1 - expcdf(y2, exponential_mean);
            elseif tag == 6
                term1 = 1 - expcdf(y1 + y2 + y3, exponential_mean);
                term2 = 1 - expcdf(y1 + y2, exponential_mean);
            elseif tag == 7
                term1 = exppdf(y3, exponential_mean);
            elseif tag == 8
                term1 = 1 - expcdf(y3, exponential_mean);
            end
        end
        if window == "W"
            if tag == 1
                term1 = exppdf(y1 + y2, exponential_mean);
                term2 = 1 - expcdf(y1, exponential_mean);
            elseif tag == 2
                term1 = exppdf(y2, exponential_mean);
            elseif tag == 3
                term1 = exppdf(y1 + y2 + y3, exponential_mean);
                term2 = 1 - expcdf(y1, exponential_mean);
            elseif tag == 4
                term1 = exppdf(y2 + y3, exponential_mean);
            elseif tag == 5
                term1 = 1 - expcdf(y2 + y3, exponential_mean);
            elseif tag == 6
                term1 = 1 - expcdf(y1 + y2 + y3, exponential_mean);
                term2 = 1 - expcdf(y1, exponential_mean);
            elseif tag == 7
                term1 = exppdf(y3, exponential_mean);
            elseif tag == 8
                term1 = 1 - expcdf(y3, exponential_mean);
            end
        end
    
        if window == "W1"
            if tag == 1
                L = term1./ term2;
            elseif tag == 2
                L = term1;
            elseif tag == 3
                L = term1./ term2;
            elseif tag == 4
                L = term1;
            elseif tag == 5
                L = term1;
            elseif tag == 6
                L = term1./ term2;
            elseif tag == 7
                L = term1;
            elseif tag == 8
                L = term1;
            end
        end
        if window == "W2"
            if tag == 1
                L = term1;
            elseif tag == 2
                L = term1;
            elseif tag == 3
                L = term1 ./ term2;
            elseif tag == 4
                L = term1 ./ term2;
            elseif tag == 5
                L = term1./ term2;
            elseif tag == 6
                L = term1./ term2;
            elseif tag == 7
                L = term1;
            elseif tag == 8
                L = term1;
            end
        end
        if window == "W"
            if tag == 1
                L = term1 ./ term2;
            elseif tag == 2
                L = term1;
            elseif tag == 3
                L = term1 ./ term2;
            elseif tag == 4
                L = term1;
            elseif tag == 5
                L = term1;
            elseif tag == 6
                L = term1./ term2;
            elseif tag == 7
                L = term1;
            elseif tag == 8
                L = term1;
            end
        end
    end   
    if ~isempty(L)
        LL = [LL; L];
    end
end

f = -sum(log( LL ));   % -log-likelihood as fminsearch is to find the minimum