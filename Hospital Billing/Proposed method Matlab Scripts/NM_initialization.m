function Theta = NM_initialization(data, m)
% data: the data to be fitting including complete, left-truncated right-censored
% window: the window to be fitted, can be W1, W2 or W
% m: number of mixing components
% h: compensation ratio for truncated and censored data
% Theta: estimated initial gamma parameters for complete, truncated and censored data

y2 = data(:,3);
y3 = data(:,4);
V = y2+y3;
    
name_cluster = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"];
% ===== starting values for Nelder-Mead by moments
idx = kmeans(V,m);
dict = struct;          % dict: a struct to store data in k clusters 
for i = [name_cluster(1:m);1:m]  % i(1) the name of cluster, i(2) the value of the cluster
    dict.(i(1)) = V(idx==str2num(i(2)));
end

mixing_proportions = [];   % mixing proportions of the k gamma
for i = name_cluster(1:m)
    mixing_proportions = [mixing_proportions; length(dict.(i))];
end
mixing_proportions = mixing_proportions / sum(mixing_proportions);

Theta = zeros(m, 3);  % parameters of gammas, mixing, shape and scale
for i = [name_cluster(1:m);1:m]
    data_cluster = dict.(i(1));
    fir = sum(data_cluster) / length(data_cluster);
    sec = sum(data_cluster .* data_cluster) / length(data_cluster);
    rate = fir / (sec - fir * fir);   % rate
    shape = rate * fir;           % shape
    Theta(str2num(i(2)),:) = [mixing_proportions(str2num(i(2))), shape, 1/rate];
end