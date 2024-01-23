function parameters = NM_optimization(data, param, Window, niters)
    % data: data to be fitted
    % param: starting values for Nelder-Mead (mixing, shape and scale (original not log)
      
%     options = optimset('Display','off');  % options for Nealder-Mead
%     options = optimset('Display','on');  % options for Nealder-Mead
%     options = optimset('Display','iter');
    options = optimset('MaxIter',niters, 'MaxFunEvals',niters,'Display','off');  % options for Nealder-Mead
   
    [m, ~] = size(param); 
    param_log = NM_param2log(param, m);              % transform param to its logrithm form
    % x is a list of parameters after optimization (logrithmized mixing, shape and mean)
    gammixfit = @(param) gammix_in(param, data, Window);
    [x, ~] = fminsearch(gammixfit, param_log, options);  
    parameters = NM_log2param(x, m);