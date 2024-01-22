function y = cdf(parameters, x) 
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
        beta = mean / alpha;     % scale

        cdf_single = gamcdf(x, alpha, beta);
        Z(i,:) = mixing * cdf_single;
    end
    y = sum(Z,1);
end