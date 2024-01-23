clear all; clc; close all;
warning off;
format long g   % not use scientific notation
% ==================================== settings ========================
Path_in = "data/sliding windows/";

nbins = 50;         % used in stage 6 for plotting histogram
fontSize = 20;
max_hours = 12*31*24*2;        % 2 years
    
ms = [2,3,4,5,6,7];   % the number of mixing components
% stage3: short-runs 
% stage4: select the best model under each sliding window with different h and m based on BIC
          % use the average likelihood of the three windows (keep all models with the same mixing components) 
% stage5: long-runs of Nelder-Mead based on models from stage 4
% stage8: likelihood ratio test
% stage10: long-runs with multi-rounds of Nelder-Mead based on models from stage 5
% stage11: cdf with Kaplan Miere from R 
% stage12: the number of complete, truncated and censored data within the three windows
% stage13: comparing W1 and W2 within the 1st and 2nd window (mean, variance)

Stages = [3,4,5,10,8,11,12];
[stage8_out, stage10_out, stage12_out] = deal([],[],[]);
stage10_rounds = 5;
startT = clock;
for stage = Stages
    if stage == 0
        % ================== transform data to 5 columns, ID, y1,y2,y3 and tag            
        months = [
         "01/01/2013 00:00"...,
         "02/01/2013 00:00"...,
         "03/01/2013 00:00"...,
         "04/01/2013 00:00"...,
         "05/01/2013 00:00"...,
         "06/01/2013 00:00"...,
         "07/01/2013 00:00"...,
         "08/01/2013 00:00"...,
         "09/01/2013 00:00"...,
         "10/01/2013 00:00"...,
         "11/01/2013 00:00"...,
         "12/01/2013 00:00"...,
         "01/01/2014 00:00"...,
         "02/01/2014 00:00"...,
         "03/01/2014 00:00"...,
         "04/01/2014 00:00"...,
         "05/01/2014 00:00"...,
         "06/01/2014 00:00"...,
         "07/01/2014 00:00"...,
         "08/01/2014 00:00"...,
         "09/01/2014 00:00"...,
         "10/01/2014 00:00"...,
         "11/01/2014 00:00"...,
         "12/01/2014 00:00"...,
         "01/01/2015 00:00"...,
         "02/01/2015 00:00"...,
         "03/01/2015 00:00"...,
         "04/01/2015 00:00"...,
         "05/01/2015 00:00"...,
         "06/01/2015 00:00"...,
         "07/01/2015 00:00"...,
         "08/01/2015 00:00"...,
         "09/01/2015 00:00"...,
         "10/01/2015 00:00"...,
         "11/01/2015 00:00"...,
         "12/01/2015 00:00"
        ];

        % windows for concept drift detection
        W1 = [];
        W2 = [];
        [~, nmonths] = size(months);
        nwindow = 12; % the length of the observation window, 12,i.e., 1 year
        for Ts1 = 1:12
            temp_w1 = [months(Ts1), months(Ts1+nwindow)];
            temp_w2 = [months(Ts1+nwindow), months(Ts1+2*nwindow)];
            W1 = [W1;temp_w1];
            W2 = [W2;temp_w2];
        end

        data_file = strcat('data/CSDM2.csv');
        data = readmatrix(data_file, 'OutputType', 'string');   % data to be fitted
        Path_output = "data/sliding windows/";
        for window = 1:12
            W1_single = W1(window,:);
            W2_single = W2(window,:);
            new_data = data_extraction(data, W1_single, W2_single);
            new_data2 = sortrows(new_data,5);

            T = array2table(new_data2);
            T.Properties.VariableNames = {'ID', 'y1', 'y2', 'y3', 'tag'};
            output_file = strcat(Path_output, num2str(window), '.csv');
            writetable(T, output_file);
        end
    end
end   
for sliding_window = 1:12
    W_file = strcat(Path_in, num2str(sliding_window), '.csv');
    data = readmatrix(W_file, 'OutputType', 'double');   % data to be fitted
    [n,~] = size(data);
    % ======================== Nelder-Mead initialization ===========
    for stage = Stages 
        [sliding_window, stage]
        if stage == 3
            n_iters = 100; % the number of iterations in Nelder-Mead
            % create folders for saving model parameters and images for all sliding windows
            Path_slidingwindow = strcat('output/models/', num2str(sliding_window));
            Path_W = strcat('output/models/', num2str(sliding_window), '/W/');    
            Path_W1 = strcat('output/models/', num2str(sliding_window), '/W1/');
            Path_W2 = strcat('output/models/', num2str(sliding_window), '/W2/');
            if exist(Path_slidingwindow, 'dir')==0 
                mkdir (Path_slidingwindow)
                mkdir ( Path_W );
                mkdir ( Path_W1 );
                mkdir ( Path_W2 );
            end
            
            for m = ms
                param_W = NM_initialization(data, m); % mixing, shape and scale
                [param_W1, param_W2] = deal(param_W, param_W);
                
                [data_good, ~, ~, ~] = NM_loglikelihood_gooddata(data, param_W, param_W1, param_W2);
                % optimize models (short runs)      
                param_shortW = NM_optimization(data_good, param_W, 'W', n_iters);
                param_shortW1 = NM_optimization(data_good, param_W1, 'W1', n_iters);
                param_shortW2 = NM_optimization(data_good, param_W2, 'W2', n_iters);
                NM_saving(param_shortW, Path_W, 'short', m);
                NM_saving(param_shortW1, Path_W1, 'short', m);
                NM_saving(param_shortW2, Path_W2, 'short', m);
            end
            
        end
        if stage == 4
            % create folders for saving best models as short run models
            Path_slidingwindow = strcat('output/short run/', num2str(sliding_window));
            if exist(Path_slidingwindow, 'dir')==0 
                mkdir (Path_slidingwindow)
            end
            
            % the path of models to be estimated
            Path_W = strcat('output/models/', num2str(sliding_window), '/W/');    
            Path_W1 = strcat('output/models/', num2str(sliding_window), '/W1/');
            Path_W2 = strcat('output/models/', num2str(sliding_window), '/W2/');  
            target_Path = strcat('output/short run/', num2str(sliding_window), '/' );    
            
            [L_temp, L1_temp, L2_temp] = deal([],[],[]);
            for m = ms
                % best model of W under different m
                W_file = strcat(Path_W, num2str(m), '.csv');
                W_file1 = strcat(Path_W1, num2str(m), '.csv');
                W_file2 = strcat(Path_W2, num2str(m), '.csv');
                parameters_W = readmatrix( W_file, 'OutputType', 'double');   % data to be fitted    
                parameters_W1 = readmatrix(W_file1, 'OutputType', 'double');   % data to be fitted
                parameters_W2 = readmatrix(W_file2, 'OutputType', 'double');   % data to be fitted

                param_logW = NM_param2log(parameters_W, m);        % transform param to its logrithm form
                param_logW1 = NM_param2log(parameters_W1, m);        % transform param to its logrithm form
                param_logW2 = NM_param2log(parameters_W2, m);        % transform param to its logrithm
                    
                [data_good, ~, ~, ~] = NM_loglikelihood_gooddata(data, parameters_W, parameters_W1, parameters_W2);
                [W_fg, ~, ~] = gammix_in(param_logW, data_good, "W");
                [W1_fg, ~, ~] = gammix_in(param_logW1, data_good, "W1");
                [W2_fg, ~, ~] = gammix_in(param_logW2, data_good, "W2");
                
                BIC_temp = log(n) * (3*m-1) - 2*W_fg;
                BIC1_temp = log(n) * (3*m-1) - 2*W1_fg;
                BIC2_temp = log(n) * (3*m-1) - 2*W2_fg; 
                
                [L_temp, L1_temp, L2_temp] = deal([L_temp, BIC_temp],[L1_temp, BIC1_temp],[L2_temp, BIC2_temp]);
            end
            [LL_W, ind_W] = max(L_temp);
            [LL_W1, ind_W1] = max(L1_temp);
            [LL_W2, ind_W2] = max(L2_temp);
               
            Best_fileW = strcat( num2str(ms(ind_W)), ".csv"); 
            Best_fileW1 = strcat( num2str(ms(ind_W1)), ".csv"); 
            Best_fileW2 = strcat( num2str(ms(ind_W2)), ".csv"); 
            best_files = [Best_fileW, Best_fileW1, Best_fileW2];
            
            % copy best models to folder best
            NM_copy(best_files, Path_W, Path_W1, Path_W2, target_Path)
        end
        if stage == 5
            n_iters = 2000; % the number of iterations in Nelder-Mead
            % create folders for saving long-run models 
            Path_longrun = strcat('output/long run/', num2str(sliding_window), '/');
            if exist(Path_longrun, 'dir')==0 
                mkdir (Path_longrun)
            end

            % load models
            parameters_Path = strcat('output/short run/', num2str(sliding_window), '/' );    
            parameters_W = readmatrix( strcat(parameters_Path, 'W.csv'), 'OutputType', 'double');                                        
            parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
            parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double');   
            
            % select good data
            [data_good, ~, ~, ~] = NM_loglikelihood_gooddata(data, parameters_W, parameters_W1, parameters_W2);
            % optimize models            
            param_longW = NM_optimization(data_good, parameters_W, 'W', n_iters);
            param_longW1 = NM_optimization(data_good, parameters_W1, 'W1', n_iters);
            param_longW2 = NM_optimization(data_good, parameters_W2, 'W2', n_iters);
            NM_saving(param_longW, Path_longrun, 'long', 'W')
            NM_saving(param_longW1, Path_longrun, 'long', 'W1')
            NM_saving(param_longW2, Path_longrun, 'long', 'W2')
        end
        if stage == 10
            n_iters = 2000;
            % create folders for saving long-run models 
            Path_longrun = strcat('output/long run2/', num2str(sliding_window), '/');
            if exist(Path_longrun, 'dir')==0 
                mkdir (Path_longrun)
            end
            for loop = 1:stage10_rounds
                % load models
                if loop == 1
                    parameters_Path = strcat('output/long run/', num2str(sliding_window), '/' );  
                else
                    parameters_Path = strcat('output/long run2/', num2str(sliding_window), '/' );  
                end
                parameters_W = readmatrix( strcat(parameters_Path, 'W.csv'), 'OutputType', 'double');                                        
                parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
                parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double');   
            
                % select good data
                [data_good, ~, ~, ~] = NM_loglikelihood_gooddata(data, parameters_W, parameters_W1, parameters_W2);
                % optimize models            
                param_longW = NM_optimization(data_good, parameters_W, 'W', n_iters);
                param_longW1 = NM_optimization(data_good, parameters_W1, 'W1', n_iters);
                param_longW2 = NM_optimization(data_good, parameters_W2, 'W2', n_iters);
                
                % remove mixing components with extremely small mixing proportions
                param_longW = NM_parameters_modification(param_longW);
                param_longW1 = NM_parameters_modification(param_longW1);
                param_longW2 = NM_parameters_modification(param_longW2);             
                NM_saving(param_longW, Path_longrun, 'long', 'W')
                NM_saving(param_longW1, Path_longrun, 'long', 'W1')
                NM_saving(param_longW2, Path_longrun, 'long', 'W2')
            end
        end
        if stage == 8
            nsamples = 140; % the number of samples for likelihood ratio test
            Path_out = 'output/test/likelihood ratio test.csv';
            % likelihood ratio test for short runs
            parameters_Path = strcat('output/short run/', num2str(sliding_window), '/' );    
            [C_short, pValue_short, dof_short] = NM_LRTest(parameters_Path, data, nsamples);

            % likelihood ratio test for long runs
            parameters_Path = strcat('output/long run/', num2str(sliding_window), '/' );    
            [C_long, pValue_long, dof_long] = NM_LRTest(parameters_Path, data, nsamples);
            
            % likelihood ratio test for long runs2
            parameters_Path = strcat('output/long run2/', num2str(sliding_window), '/' );    
            [C_long2, pValue_long2, dof_long2] = NM_LRTest(parameters_Path, data, nsamples);

            stage8_out = [stage8_out; [sliding_window, C_short, pValue_short, dof_short, C_long, ...
                pValue_long, dof_long, C_long2, pValue_long2, dof_long2]];
            T = array2table(stage8_out);
%             T.Properties.VariableNames(1:3) = {'sliding window', 'short', 'long'};
            writetable(T, Path_out);
        end
        if stage == 11
            % cdf and KM for short runs
            path_R = strcat( 'data/KMR/', num2str(sliding_window) );
            parameters_Path = strcat('output/short run/', num2str(sliding_window), '/' );    
            output_png = strcat('output/CDF2/short run/', num2str(sliding_window), '.png' );
            NM_KM_R(parameters_Path, data, max_hours, output_png, path_R);
            
            % cdf and KM for long runs
            parameters_Path = strcat('output/long run/', num2str(sliding_window), '/' );    
            output_png = strcat('output/CDF2/long run/', num2str(sliding_window), '.png' );
            NM_KM_R(parameters_Path, data, max_hours, output_png, path_R)

            % cdf and KM for long runs2
            parameters_Path = strcat('output/long run2/', num2str(sliding_window), '/' );    
            output_png = strcat('output/CDF2/long run2/', num2str(sliding_window), '.png' );
            NM_KM_R(parameters_Path, data, max_hours, output_png, path_R)
        end
        if stage == 12
            Path_out = 'output/test/summary.csv';
            tags = data(:,5);
            % =================== window W
            tags_W_complete = tags(tags == 2 |tags == 4 |tags == 7);
            tags_W_LT = tags(tags == 1 |tags == 3);
            tags_W_RC = tags(tags == 5 |tags == 8);
            tags_W_LTRC = tags(tags == 6);
            [n_W_complete, ~] = size(tags_W_complete); 
            [n_W_LT, ~] = size(tags_W_LT); 
            [n_W_RC, ~] = size(tags_W_RC); 
            [n_W_LTRC, ~] = size(tags_W_LTRC);  

            % ===================== window W1            
            tags_W1_complete = tags(tags == 2);
            tags_W1_LT = tags(tags == 1);
            tags_W1_RC = tags(tags == 4 |tags == 5);
            tags_W1_LTRC = tags(tags == 3 | tags == 6);
            [n_W1_complete, ~] = size(tags_W1_complete); 
            [n_W1_LT, ~] = size(tags_W1_LT); 
            [n_W1_RC, ~] = size(tags_W1_RC); 
            [n_W1_LTRC, ~] = size(tags_W1_LTRC); 
            
            % ===================== window W2           
            tags_W2_complete = tags(tags == 7);
            tags_W2_LT = tags(tags == 3 |tags == 4);
            tags_W2_RC = tags(tags == 8);
            tags_W2_LTRC = tags(tags == 5 | tags == 6);
            [n_W2_complete, ~] = size(tags_W2_complete); 
            [n_W2_LT, ~] = size(tags_W2_LT); 
            [n_W2_RC, ~] = size(tags_W2_RC); 
            [n_W2_LTRC, ~] = size(tags_W2_LTRC);  
            
            n_W1 = n_W1_complete + n_W1_LT + n_W1_RC + n_W1_LTRC;
            n_W2 = n_W2_complete + n_W2_LT + n_W2_RC + n_W2_LTRC;       
            stage12_out = [stage12_out; [sliding_window, ...
                round(n_W_complete/n,3), round(n_W_LT/n,3), round(n_W_RC/n,3), round(n_W_LTRC/n,3), ...
                round(n_W1_complete/n_W1,3), round(n_W1_LT/n_W1,3), round(n_W1_RC/n_W1,3), round(n_W1_LTRC/n_W1,3), ...
                round(n_W2_complete/n_W2,3), round(n_W2_LT/n_W2,3), round(n_W2_RC/n_W2,3), round(n_W2_LTRC/n_W2,3)]];
            T = array2table(stage12_out);
            writetable(T, Path_out);            
        end
        if stage == 13
            parameters_Path = strcat('output/long run2/', num2str(sliding_window), '/' );  
            parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
            parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double'); 
            parameters_W1(:,4) = parameters_W1(:,2) .* parameters_W1(:,3);   % calculate means
            parameters_W2(:,4) = parameters_W2(:,2) .* parameters_W2(:,3);   % calculate means
            
%             parameters_W1 = sortrows(parameters_W1,-1);
%             parameters_W2 = sortrows(parameters_W2,-1);
            
            data_al  = data(:,2) + data(:,3) + data(:,4);
            samples = 1:100:max_hours;
            y1 = pdf(parameters_W1(:,[1,2,4]), samples);
            y2 = pdf(parameters_W2(:,[1,2,4]), samples);
            
            fig = figure;
            nbins = 25;
            histogram(data_al, nbins, 'Normalization','pdf');
            hold on;
            plot(samples, y1, 'r-', 'DisplayName', 'KM W1');
            hold on;
            plot(samples, y2, 'g-', 'DisplayName', 'KM W2');
            legend('Location', 'Best');
            output_png = strcat('output/PDF/', num2str(sliding_window), '.png' );
            print(fig,'-r500','-dpng',output_png);
            close(fig)
                
            out1 = [parameters_W1(:,1), parameters_W1(:,4) / 24]
            out2 = [parameters_W2(:,1), parameters_W2(:,4) / 24]
        end
    end
end

% data_file = strcat('output/test/summary.csv');
% data = readmatrix(data_file, 'OutputType', 'double');   % data to be fitted
% avg_summ = mean(data);
endT = clock;
Texcution = etime(endT,startT);