clear all; clc; close all;
warning off;
format long g   % not use scientific notation
% ==================================== settings ========================

% stage1: test1-type1 error on complete data
% stage2: test2-generate data
% stage3: test2-transform ID, stime, etime to ID, y1, y2,y3,tag
% stage4: test2-data fitting to exponential distributions
% stage5: test2-LR test
% stage6: test2-type1 error on left truncated and right censored data
% stage7: test2-the proportion of complete data within each window

% stages = [1,2,3,4,5];
% stages = [1];
% stages = [2,3,4,5,6,7];
stages = [1,2,3,4,5,6,7];
% ==================================== method validation ========================
startT = clock;

for stage = stages
    windows = [ 2,  5,  7, 10, 15, 20, 30, 40, 50, 60, 70, 80];
    n_run = 1000;
    n_max = 50;  % the number of time points, also the range of the histogram
    timepoints = 0:0.5:n_max;        % time points to generate n_sample samples each
    if stage == 1
        % ============================= settings
        n_sample = 1000;
        n_run = 1000;
        lambdas = 0.1:0.1:1;
        alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
        dof = 1;                     % the degree freedom
        Path_in = 'test1/';

        % ============================= LR test
        for lam = lambdas
            test1_output = [];
            for run = 1:n_run
                % generate source exponential data with mean=1/lam
                sample1 = exprnd(1 / lam, n_sample, 1);
                sample2 = exprnd(1 / lam, n_sample, 1);
                data = [sample1; sample2];

                % moments to calculate lam for W, W1 and W2
                [lamW1, lamW2] = deal( 1 / mean(sample1), 1 / mean(sample2) );
                lamW = 1 / mean(data);        

                % likelihoods for W, W1 and W2
                LC_s = exppdf(data, 1 / lamW);         % simpler model 
                LC_W1 = exppdf(sample1, 1 / lamW1);
                LC_W2 = exppdf(sample2, 1 / lamW2);
                LC_c = [LC_W1; LC_W2];  % complex model

                L_c = sum(log( LC_c ));
                L_s = sum(log( LC_s ));
                U_LL = max([L_c, L_s]);
                L_LL = min([L_c, L_s]);            
                [~,pValue] = lratiotest(U_LL, L_LL, dof);

                test1_output = [test1_output; [lamW, lamW1, lamW2, U_LL, L_LL, pValue]];
            end
            Path_out = strcat( Path_in, num2str(lam), '.csv'); 
            T = array2table(test1_output);
            writetable(T, Path_out);
        end
        % =============================  type1 error
        output = [];
        for lam = lambdas
            W_file = strcat( Path_in, num2str(lam), '.csv' );
            data = readmatrix( W_file, 'OutputType', 'double');
            pvalue = data(:,end);

            temp = [];
            for alpha = alphas
                type1 = sum( pvalue <= alpha ) / length(pvalue);
                temp = [temp, type1];
            end   
            output = [output; [lam, temp] ];
        end            
        Path_out = strcat( Path_in, 'summary.csv'); 
        T = array2table(output);
        writetable(T, Path_out);
    end
    if stage == 2
        % ============================= settings
        % fix the lambda, given different censoring time to simulate different percentage of censoring
        n_sample = 10;
        Path_in = 'test2/data/';
        lam = 0.2;
        exp_mean = 1 / lam;

        % ============================= generate source data
        for run = 1:n_run
            out = [];         % the start and end time of the samples
            for i = timepoints
                sample = exprnd(exp_mean, n_sample, 1); 
                temp = [ones(n_sample,1)*i, ones(n_sample,1)*i + sample];
                out = [out; temp];
            end
            [n_samples, ~] = size(out);
            out = [ [1:n_samples]', out ];
            Path_out = strcat( Path_in, num2str(run), '.csv'); 
            T = array2table(out);
            T.Properties.VariableNames = {'ID', 'stime', 'etime'};
            writetable(T, Path_out);            
        end
    end
    if stage == 3
        % ============================= settings
        % fix the lambda, given different censoring time to simulate different percentage of censoring
        Path_in = 'test2/data/';
        Path_out = "test2/data2/";
        % ============================= generate source data
        for window = windows
            W1 = [n_max / 2 - window / 2, n_max / 2];  % boundary of W1 
            W2 = [n_max / 2, n_max / 2 + window / 2];  % boundary of W2

            for run = 1:n_run
                data = readmatrix( strcat(Path_in, num2str(run),'.csv'), 'OutputType', 'double');  
                IDs = unique(data(:,1));  

                new_data = data_extraction(data, W1, W2);
                new_data2 = sortrows(new_data,5);

                T = array2table(new_data2);
                T.Properties.VariableNames = {'ID', 'y1', 'y2', 'y3', 'tag'};
                output_file = strcat(Path_out, num2str(window), '_', num2str(run), '.csv');
                writetable(T, output_file);  
            end
        end
    end
    if stage == 4
        % ============================= settings
        % Fitting exponential models
        Path_in = 'test2/data2/';
        Path_out = "test2/models/";
        % ============================= generate source data
        for window = windows
            parameters_out = [];   % saving lambda for W, W1 and W2 under different run
            for run = 1:n_run
                file_name = strcat(Path_in, num2str(window), '_', num2str(run), '.csv');
                data = readmatrix( file_name, 'OutputType', 'double');  
                lam_W_ini = NM_initialization(data); 
                [lam_W1_ini, lam_W2_ini] = deal(lam_W_ini, lam_W_ini);

                data_good = NM_loglikelihood_gooddata(data, lam_W_ini, lam_W1_ini, lam_W2_ini);
                % optimize models            
                out_lam = NM_optimization(data_good, lam_W_ini, 'W');
                out_lam1 = NM_optimization(data_good, lam_W1_ini, 'W1');
                out_lam2 = NM_optimization(data_good, lam_W2_ini, 'W2');

                parameters_out = [parameters_out; [run, out_lam, out_lam1, out_lam2]];           
            end
            T = array2table(parameters_out);
            T.Properties.VariableNames = {'run', 'lam', 'lam1', 'lam2'};
            output_file = strcat(Path_out, num2str(window), '.csv');
            writetable(T, output_file);
        end
    end
    if stage == 5
        % ============================= settings
        % Fitting exponential models
        dof = 1;                     % the degree freedom
        Path_out = 'test2/';
        Path_data = 'test2/data2/';
        Path_model = 'test2/models/';
        % ============================= LR test
        for window = windows
            test2_output = [];
            for run = 1:n_run
                file_name = strcat(Path_data, num2str(window), '_', num2str(run), '.csv');
                data = readmatrix( file_name, 'OutputType', 'double');  

                file_models = strcat(Path_model, num2str(window), '.csv');
                parameters = readmatrix( file_models, 'OutputType', 'double');

                [lam, lam1, lam2] = deal( parameters(run,2),parameters(run,3),parameters(run,4) );
                data_good = NM_loglikelihood_gooddata(data, lam, lam1, lam2);
                % likelihoods for W, W1 and W2
                [L_s, LL] = exp_in(log(lam), data_good, 'W');          % simpler model 
                [L_W1, LL1] = exp_in(log(lam1), data_good, 'W1');
                [L_W2, LL2] = exp_in(log(lam2), data_good, 'W2');
                L_c = L_W1 + L_W2;                  % complex model

                [n_LL,~] = size(LL);
                U_LL = max([L_c, L_s]) / n_LL * 1000;
                L_LL = min([L_c, L_s]) / n_LL * 1000;            
                [~,pValue] = lratiotest(U_LL, L_LL, dof);

                test2_output = [test2_output; [lam, lam1, lam2, U_LL, L_LL, pValue]];          
            end
            T = array2table(test2_output);
            T.Properties.VariableNames = {'lam', 'lam1', 'lam2', 'U', 'L', 'pvalue'};
            output_file = strcat(Path_out, num2str(window), '.csv');
            writetable(T, output_file);
        end
    end
    if stage == 6
        % ============================= settings
        % Fitting exponential models
        alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];                  % the degree freedom
        Path_out = 'test2/';
        Path_summary = strcat( Path_out, 'summary.csv'); 
        % =============================  type1 error
        output = [];
        for window = windows
            W_file = strcat( Path_out, num2str(window), '.csv' );
            data = readmatrix( W_file, 'OutputType', 'double');
            pvalue = data(:,end);

            temp = [];
            for alpha = alphas
                type1 = sum( pvalue <= alpha ) / length(pvalue);
                temp = [temp, type1];
            end   
            output = [output; [window, temp] ];
        end            
        T = array2table(output);
        writetable(T, Path_summary);    
    end
    if stage == 7
        % ============================= settings
        % Fitting exponential models
        Path_in = 'test2/data2/';
        Path_out = 'test2/summary2.csv';
        % =============================  type1 error
        output = [];
        for window = windows
            temp_window = [];
            for run = 1:n_run
                W_file = strcat( Path_in, num2str(window), '_', num2str(run), '.csv' );
                data = readmatrix( W_file, 'OutputType', 'double');
                [n, ~] = size(data);
                [tags_W, tags_W1, tags_W2] = deal( data(:,5),data(:,5),data(:,5) );

                tags_W(tags_W==1 | tags_W==3 | tags_W==5 | tags_W==6 | tags_W==8) = 0;
                tags_W(tags_W==2 | tags_W==4 | tags_W==7) = 1;
                tags_W1(tags_W1==1) = 0;
                tags_W1(tags_W1==2) = 1;
                tags_W2(tags_W2==1) = 0;
                tags_W2(tags_W2==7) = 1;
                
                n_complete_W = sum(tags_W == 1) / n;
                n_complete_W1 = sum(tags_W1 == 1) / n;
                n_complete_W2 = sum(tags_W2 == 1) / n;
                temp_window = [temp_window; [n_complete_W, n_complete_W1, n_complete_W2]];
            end
            output = [output; [window, mean(temp_window)]];
            T = array2table(output);
            writetable(T, Path_out);   
        end 
    end
end
endT = clock;
Texcution = etime(endT,startT);