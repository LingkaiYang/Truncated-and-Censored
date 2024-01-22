function param_log = NM_param2log(param, m)
    % param: real parmaters in matrix m*3, mixing, shape, scale
    % m: the number of gammas to mix
    % ===================== transform param to its logrithm form ==========
    parameters = param;                                     % mixing, shape, scale
    parameters(:,4) = parameters(:,2) .* parameters(:,3);   % calculate means
    parameters(:,5) = parameters(:,2) ./ parameters(:,4);  % calculate the rate of each gamma (beta)
    parameters = parameters(:,[1,2,4]); % m*3, for mixing, shape and mean 
 
    new_parameters = parameters;
    new_parameters(:,2) = log(parameters(:,2));    
    new_parameters(:,3) = log(parameters(:,3));    
    p = 1;
    for i=1:(m-1)        % mixing proportion for each gamma
        temp = parameters(i,1);
        r = log(temp/(p-temp));
        p=p-temp;
        new_parameters(i,1) = r;
    end
    new_parameters(m,1) = p;
    new_parameters = new_parameters.';      % transpose as 3*m

    param_log = new_parameters(:);   % transform as a list, 1,2,3 the first gamma,...
    param_log(3*m-2) = []; % remove the last mixing proportion to have 3m-1 free parameters

