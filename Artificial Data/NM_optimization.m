function parameters = NM_optimization(data, lam, Window)
    % data: data to be fitted
    % lam: starting values for Nelder-Mead (original not log)
      
%     options = optimset('Display','off');  % options for Nealder-Mead
%     options = optimset('Display','on');  % options for Nealder-Mead
%     options = optimset('Display','iter');
    niters = 10000;
    options = optimset('MaxIter',niters, 'MaxFunEvals',niters,'Display','off');  % options for Nealder-Mead
   
    lam_log = log(lam);              % transform param to its logrithm form
    % x is a list of parameters after optimization (logrithmized mixing, shape and mean)
    gammixfit = @(param) exp_in(param, data, Window);
    [x, ~] = fminsearch(gammixfit, lam_log, options);  
    parameters = exp(x);