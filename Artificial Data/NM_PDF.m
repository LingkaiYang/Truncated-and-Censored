function NM_PDF(parameters_Path, data, max_hours, output_png)   
    fontSize = 20;
    LineWidth = 2;
    nbins = 50;    
    parameters_W = readmatrix( strcat(parameters_Path, 'W.csv'), 'OutputType', 'double');                                  
    parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
    parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double');   
    [~, Ls, Lc1, Lc2] = NM_loglikelihood_gooddata(data, parameters_W, parameters_W1, parameters_W2);

    parameters_W(:,4) = parameters_W(:,2) .* parameters_W(:,3);   % calculate means
    parameters_W1(:,4) = parameters_W1(:,2) .* parameters_W1(:,3);   % calculate means
    parameters_W2(:,4) = parameters_W2(:,2) .* parameters_W2(:,3);   % calculate means
    % plot and save as png
    d = data(:,2) + data(:,3);
    samples = linspace(1, max_hours, 5000);
    y_W = pdf(parameters_W(:,[1,2,4]), samples);
    y_W1 = pdf(parameters_W1(:,[1,2,4]), samples);
    y_W2 = pdf(parameters_W2(:,[1,2,4]), samples);

    fig = figure;
    histogram(d, nbins,'Normalization','pdf');  
    hold on;
    plot(samples,y_W, 'DisplayName', strcat('W: ', num2str(Ls)), 'Linewidth', LineWidth)
    hold on;
    plot(samples,y_W1, 'DisplayName', strcat('W1: ', num2str(Lc1)), 'Linewidth', LineWidth)
    hold on;
    plot(samples,y_W2, 'DisplayName', strcat('W2: ', num2str(Lc2)), 'Linewidth', LineWidth)

    set(gca,'FontSize',fontSize);  % fontsize for axes and legend
    set(gca,'yticklabel',[])     % not show axis
    set(gca,'ytick',[])          % not show values
    box off
    legend('Location', 'Best');
    xlabel('Hours','FontSize',fontSize)
    ylabel('pdf','FontSize',fontSize)     

    set(fig,'visible','off');
    set(gca,'LooseInset',get(gca,'TightInset'));
    print(fig,'-r500','-dpng',output_png);
    close(fig) 
end
               
    