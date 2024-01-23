function new_parameters = NM_parameters_modification(parameters)
    % parameters: mixing, shape and scale
    [m,~] = size(parameters);
    eps = 10e-4;   % if mixing proportions are smaller than eps, the gamma is filtered
    new_parameters = parameters;
    
    ID_delete = [];
    for i = 1:m
        [mixing, shape, scale] = deal(parameters(i,1), parameters(i,2), parameters(i,3));
        if mixing < eps || isinf(shape) || isnan(shape) || isinf(scale) || isnan(scale) || isinf(mixing) || isnan(mixing)
            ID_delete = [ID_delete,i];
        end
    end
    new_parameters(ID_delete,:) = [];
    
    mixing_sum = sum(new_parameters(:,1));
	new_parameters(:,1) = new_parameters(:,1) / mixing_sum;
    
    