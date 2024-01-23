function NM_saving(parameters, Path_out, tag, args)
    % parameters: mixing, shape and scale
    T = array2table(parameters);
    T.Properties.VariableNames(1:3) = {'weight', 'shape', 'scale'};
    
    if tag == "short"
        m = args;
        output_file = strcat(Path_out, num2str(m), '.csv');
        writetable(T, output_file);
    end
    
    if tag == "long"
        window = args;
        output_file = strcat(Path_out, window, '.csv');
        writetable(T, output_file);
    end

    
    