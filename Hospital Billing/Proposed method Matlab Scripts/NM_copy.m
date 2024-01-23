function NM_copy(best_files, Path_W, Path_W1, Path_W2, target_Path)
    % copy best models to folder best
    bestW = best_files(1);
    bestW1 = best_files(2);
    bestW2 = best_files(3);
    
    sourceW = char( fullfile(Path_W, bestW) );
    targetW = strcat(target_Path, 'W.csv');
    copyfile(sourceW, targetW)  

    sourceW1 = char( fullfile(Path_W1, bestW1) );
    targetW1 = strcat(target_Path, 'W1.csv');
    copyfile(sourceW1, targetW1)  
    
    sourceW2 = char( fullfile(Path_W2, bestW2) );
    targetW2 = strcat(target_Path, 'W2.csv');
    copyfile(sourceW2, targetW2)  
end
               
    