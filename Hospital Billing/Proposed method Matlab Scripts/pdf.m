function y = pdf(parameters, x) 
    % parameters: gamma mix class with size m*3, mixing, shape and mean
    % x: data to calculate the pdf 
%     k = length(parameters);    % error when a single gamma, use size() instead
    [k,~] = size(parameters);
    n = length(x);
    Z = zeros(k,n);    % logpdf
    for i=1:k        % log(pdf) for each gamma component
        mixing = parameters(i,1);
        alpha = parameters(i,2);
        mean = parameters(i,3);
        beta = alpha / mean;

        logpdf = alpha * log(beta) + (alpha-1)*log(x)-beta*x-gammaln(alpha);
        Z(i,:) = mixing * exp(logpdf);
    end
    y = sum(Z,1);
end