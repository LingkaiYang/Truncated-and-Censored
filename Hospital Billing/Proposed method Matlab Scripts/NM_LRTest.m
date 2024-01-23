function [C, pValue, dof] = NM_LRTest(parameters_Path, data, nsamples)
        [n, ~] = size(data);
        % likelihood ratio test
%         parameters_Path = strcat('output/long run2/', num2str(sliding_window), '/' );    
        parameters_W = readmatrix( strcat(parameters_Path, 'W.csv'), 'OutputType', 'double');                                        % mixing, shape, scale
        parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
        parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double');           
        [m, ~] = size(parameters_W);
        [m1, ~] = size(parameters_W1);
        [m2, ~] = size(parameters_W2);

        [~, Ls, Lc1, Lc2] = NM_loglikelihood_gooddata(data, parameters_W, parameters_W1, parameters_W2);
        Lc = Lc1+Lc2;
        C = abs(Lc-Ls);
        dof = (m1+m2-m)*3-1;
%         dof = m*3-1;

        Lc_shrink_long = Lc / n * nsamples;
        Ls_shrink_long = Ls / n * nsamples;
        U_LL = max([-Lc_shrink_long,-Ls_shrink_long]);
        L_LL = min([-Lc_shrink_long,-Ls_shrink_long]);
        [~,pValue] = lratiotest(U_LL, L_LL, dof);