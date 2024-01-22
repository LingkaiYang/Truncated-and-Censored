function data_good = NM_loglikelihood_gooddata(data, lam, lam1, lam2)
    % select good data for log-likelihood calculation
    IDSS = data(:,1);
    lam_log = log(lam);              % transform param to its logrithm form
    lam1_log = log(lam1);              % transform param to its logrithm form
    lam2_log = log(lam2);              % transform param to its logrithm form

    [~, W_L] = exp_in(lam_log, data, "W");
    [~, W1_L] = exp_in(lam1_log, data, "W1");
    [~, W2_L] = exp_in(lam2_log, data, "W2");

    ind_W = NM_Nan_Inf_0(IDSS, W_L);
    ind_W1 = NM_Nan_Inf_0(IDSS, W1_L);
    ind_W2 = NM_Nan_Inf_0(IDSS, W2_L);
    ID_bad = [ind_W, ind_W1, ind_W2]; 

    IDs_good = setdiff(IDSS,ID_bad).';    % case IDs of good instances
    ind_good = [];                        % indices of good instances
    for id = IDs_good
        temp = find(IDSS==id);
        ind_good = [ind_good, temp];
    end

    data_good = data(ind_good,:); 
end
               
    