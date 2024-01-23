function NM_KM_R(parameters_Path, data, max_hours, output_png, Path_R)
    % ===================================== cdf and kaplan-miere
    fontSize = 20;
    LineWidth = 2;
    % cdf of instances under W, W1 and W2    
    parameters_W = readmatrix( strcat(parameters_Path, 'W.csv'), 'OutputType', 'double');     % mixing, shape, scale
    parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
    parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double');   
    [~, Ls, Lc1, Lc2] = NM_loglikelihood_gooddata(data, parameters_W, parameters_W1, parameters_W2);

    parameters_W(:,4) = parameters_W(:,2) .* parameters_W(:,3);   % calculate means
    parameters_W1(:,4) = parameters_W1(:,2) .* parameters_W1(:,3);   % calculate means
    parameters_W2(:,4) = parameters_W2(:,2) .* parameters_W2(:,3);   % calculate means            

    samples = linspace(1, max_hours, max_hours);
    y_W = cdf(parameters_W(:,[1,2,4]), samples);       % cdf of samples under model W
    y_W1 = cdf(parameters_W1(:,[1,2,4]), samples);
    y_W2 = cdf(parameters_W2(:,[1,2,4]), samples);    

    % kaplan-miere estimation for cdf based on W
    km_data = readmatrix( strcat(Path_R, '_W.csv'), 'OutputType', 'double');
    km_x = km_data(:,1);
    km_cdf = km_data(:,2);
    km_cdf(end+1) = km_cdf(end);
    km_x(end+1) = samples(end);

    % kaplan-miere estimation for cdf based on W1
    km_data1 = readmatrix( strcat(Path_R, '_W1.csv'), 'OutputType', 'double');
    km_x1 = km_data1(:,1);
    km_cdf1 = km_data1(:,2);
    km_cdf1(end+1) = km_cdf1(end);
    km_x1(end+1) = samples(end);

    % kaplan-miere estimation for cdf based on W2
    km_data2 = readmatrix( strcat(Path_R, '_W2.csv'), 'OutputType', 'double');
    km_x2 = km_data2(:,1);
    km_cdf2 = km_data2(:,2);
    km_cdf2(end+1) = km_cdf2(end);
    km_x2(end+1) = samples(end);

    fig = figure;
    plot(km_x, km_cdf, 'b-', 'DisplayName', 'KM W', 'Linewidth', LineWidth)
    hold on;
    plot(km_x1, km_cdf1, 'g-', 'DisplayName', 'KM W1', 'Linewidth', LineWidth)
    hold on;
    plot(km_x2, km_cdf2, 'r-', 'DisplayName', 'KM W2', 'Linewidth', LineWidth)
    hold on;
    plot(samples,y_W, 'b--', 'DisplayName', strcat('W: ', num2str(Ls)), 'Linewidth', LineWidth)
    hold on;
    plot(samples,y_W1, 'g--', 'DisplayName',strcat('W1: ', num2str(Lc1)), 'Linewidth', LineWidth)
    hold on;
    plot(samples,y_W2, 'r--', 'DisplayName',strcat('W2: ', num2str(Lc2)), 'Linewidth', LineWidth)

    set(gca,'FontSize',fontSize);  % fontsize for axes and legend
    set(gca,'yticklabel',[])     % not show axis
    set(gca,'ytick',[])          % not show values
    box off
    legend('Location', 'Best');
    xlabel('Days','FontSize',fontSize)
    ylabel('cdf','FontSize',fontSize)    

    xticks([0 max_hours*0.25 max_hours*0.5 max_hours*0.75 max_hours])
    max_days = max_hours /24;
    xticklabels({0, int32(max_days*0.25), int32(max_days*0.5), int32(max_days*0.75), int32(max_days)})
    set(fig,'visible','off');   

    set(gca,'LooseInset',get(gca,'TightInset'));
    print(fig,'-r500','-dpng',output_png);
    close(fig)     
end
               
    