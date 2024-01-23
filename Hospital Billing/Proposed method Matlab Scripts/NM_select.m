function Best_fileW = NM_select(Path_W, m, ind_W)
    % extract all csv files for W/W1/W2 and copy the best model to folder best
    D_W = dir( Path_W );   % all csv files under m components
    files_W = cell( length(D_W),1 );   
    for i = 1:length(D_W)
        files_W(i) = {D_W(i).name};
    end

    Best_fileW = files_W{ind_W};
%     sourceW = fullfile(Path_W, Best_fileW);
%     targetW = strcat(target_Path, Window, '.csv');
%     copyfile(sourceW, targetW)    
    