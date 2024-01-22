function lam = NM_initialization(data)
% initialize W, W1 and W2 using y2+y3 based on moments
% data: the data to be fitting including complete, left-truncated right-censored
% window: the window to be fitted, can be W1, W2 or W
% Theta: estimated initial parameter

y2 = data(:,3);
y3 = data(:,4);
d = y2+y3;
    
lam = 1 / mean(d);

end