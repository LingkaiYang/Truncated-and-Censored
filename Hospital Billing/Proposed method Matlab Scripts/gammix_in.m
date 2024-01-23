function [f, Pi,LL] = gammix_in(param, data, window)
% Likelihood contribution of mixtures of gammas for incomplete 
% (letf-truncated and right-censored data), 
% parameters are a list of logrized mixing, shape and mean

% data: the data to be fitting with 5 columns, ID, y1,y2,y3 and tag
% param: gamma parameters with m mcomponents for mixing, shape and mean
%        as a list with m shapes and means and m-1 mixing components
% window: the window to be fitted, can be W1, W2 or W

%================================ preprocess initial parameters ============  
m = int8(length(param) / 3);   % the number of components
Pi = zeros(m,1);   % mixing proportions

p = 1;
for i=1:(m-1)        % mixing proportion for each gamma
    temp = exp(param(3*(i-1)+1));
    pi = p * temp / (1+temp);
    p = p-pi;
    Pi(i) = pi;
end
Pi(m) = p;

LL = [];   % loglikelihood of each gamma
for tag=1:8
    ind = data(:,5) == tag;  % indices to elements in 5th column of data that satisfy the equality
    temp_data = data(ind,:);
    [temp_n, ~] = size(temp_data); 
    
    term1 = zeros(m,temp_n);
    term2 = zeros(m,temp_n);
    for i=1:m        % log(pdf) for each gamma component
        % transform the list parameters back to shape and rate (alpha, beta)
        if i ~= m
            alpha = exp(param(3*(i-1)+2));      % shape parameter
            mean = exp(param(3*(i-1)+3));
            beta = mean / alpha;              % scale parameter
        else
            alpha = exp(param(3*(i-1)+1));
            mean = exp(param(3*(i-1)+2)); 
            beta = mean / alpha;
        end
        % calculate pdf and survival function according to the 8 scenarios
        if temp_n ~= 0
            y1 = temp_data(:,2);
            y2 = temp_data(:,3);
            y3 = temp_data(:,4);
            
%             if i == 1
%                             [y1,y2,y3]
%             end
            if window == "W1"
                if tag == 1
                    term1(i,:) = gampdf(y1 + y2, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1, alpha, beta);
                elseif tag == 2
                    term1(i,:) = gampdf(y2, alpha, beta);
                elseif tag == 3
                    term1(i,:) = 1 - gamcdf(y1 + y2, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1, alpha, beta);
                elseif tag == 4
                    term1(i,:) = 1 - gamcdf(y2, alpha, beta);
%                     [window, S]
                elseif tag == 5
                    term1(i,:) = 1 - gamcdf(y2, alpha, beta);
                elseif tag == 6
                    term1(i,:) = 1 - gamcdf(y1 + y2, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1, alpha, beta);
                elseif tag == 7
                    term1(i,:) = ones(temp_n,1);
                elseif tag == 8
                    term1(i,:) = ones(temp_n,1);
                end
            end
            if window == "W2"
                if tag == 1
                    term1(i,:) = ones(temp_n,1);
                elseif tag == 2
                    term1(i,:) = ones(temp_n,1);
                elseif tag == 3
                    term1(i,:) = gampdf(y1 + y2 + y3, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1 + y2, alpha, beta);
                elseif tag == 4
                    term1(i,:) = gampdf(y2 + y3, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y2, alpha, beta);
%                     [window, pdf, S, L]
                elseif tag == 5
                    term1(i,:) = 1 - gamcdf(y2 + y3, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y2, alpha, beta);
                elseif tag == 6
                    term1(i,:) = 1 - gamcdf(y1 + y2 + y3, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1 + y2, alpha, beta);
                elseif tag == 7
                    term1(i,:) = gampdf(y3, alpha, beta);
                elseif tag == 8
                    term1(i,:) = 1 - gamcdf(y3, alpha, beta);
                end
            end
            if window == "W"
                if tag == 1
                    term1(i,:) = gampdf(y1 + y2, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1, alpha, beta);
                elseif tag == 2
                    term1(i,:) = gampdf(y2, alpha, beta);
                elseif tag == 3
                    term1(i,:) = gampdf(y1 + y2 + y3, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1, alpha, beta);
                elseif tag == 4
                    term1(i,:) = gampdf(y2 + y3, alpha, beta);
                elseif tag == 5
                    term1(i,:) = 1 - gamcdf(y2 + y3, alpha, beta);
                elseif tag == 6
                    term1(i,:) = 1 - gamcdf(y1 + y2 + y3, alpha, beta);
                    term2(i,:) = 1 - gamcdf(y1, alpha, beta);
                elseif tag == 7
                    term1(i,:) = gampdf(y3, alpha, beta);
                elseif tag == 8
                    term1(i,:) = 1 - gamcdf(y3, alpha, beta);
                end
            end
        end
    end
    
    for i=1:m        % 
        term1(i,:) = Pi(i)*term1(i,:);
        term2(i,:) = Pi(i)*term2(i,:);
    end
    
    term1_f = sum(term1,1);
    term2_f = sum(term2,1);
    if window == "W1"
        if tag == 1
            L = term1_f./ term2_f;
        elseif tag == 2
            L = term1_f;
        elseif tag == 3
            L = term1_f./ term2_f;
        elseif tag == 4
            L = term1_f;
%                     [window, S]
        elseif tag == 5
            L = term1_f;
        elseif tag == 6
            L = term1_f./ term2_f;
        elseif tag == 7
            L = term1_f;
        elseif tag == 8
            L = term1_f;
        end
    end
    if window == "W2"
        if tag == 1
            L = term1_f;
        elseif tag == 2
            L = term1_f;
        elseif tag == 3
            L = term1_f ./ term2_f;
        elseif tag == 4
            L = term1_f ./ term2_f;
        elseif tag == 5
            L = term1_f./ term2_f;
        elseif tag == 6
            L = term1_f./ term2_f;
        elseif tag == 7
            L = term1_f;
        elseif tag == 8
            L = term1_f;
        end
    end
    if window == "W"
        if tag == 1
            L = term1_f ./ term2_f;
        elseif tag == 2
            L = term1_f;
        elseif tag == 3
            L = term1_f ./ term2_f;
        elseif tag == 4
            L = term1_f;
        elseif tag == 5
            L = term1_f;
        elseif tag == 6
            L = term1_f./ term2_f;
        elseif tag == 7
            L = term1_f;
        elseif tag == 8
            L = term1_f;
        end
    end
    if ~isempty(L)
        LL = [LL, L];
    end
end

f = -sum(log( LL ));   % -log-likelihood as fminsearch is to find the minimum