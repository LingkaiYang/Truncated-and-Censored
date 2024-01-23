function parameters = NM_log2param(x, m)
    % x: parameters after optimization, i.e., a list of logrithmized mixing, shape and mean
    % m: the number of gammas to mix
    
    % parameters: parameters in matrix m*3 for mixing, shape and scale
    % ========================= change the optimised parameters back to gamma,
    list_x = exp(x);
    p = 1;
    for i=1:(m-1)        % mixing proportion for each gamma
        temp = list_x(3*(i-1)+1);
        pi = p * temp / (1+temp);
        list_x(3*(i-1)+1) = pi;
        p = p-pi;
    end
    list_x = [list_x(1:(m-1)*3);p;list_x((m-1)*3+1:m*3-1)];

    % ======================================= change the list back to the matrix
    matrix_x = reshape(list_x,3,m).';   % columns for mixing, shape and mean

    % ======================================= calculate the scale parmeters
    scales = matrix_x(:,3) ./ matrix_x(:,2);
    matrix_x(:,4) = scales;
    parameters = matrix_x(:,[1,2,4]);

