function new_data = data_extraction(data, W1, W2) 
    % data: original dataset, ID, start time and end time
    % W1: the start and end of the 1st window
    % W2: the start and end of the 2nd window
    Ts1_time = W1(1);
    Te1_time = W1(2);
    Ts2_time = W2(1);
    Te2_time = W2(2);
    
    [n,~] = size(data);
    new_data = zeros(n,5);

    for i = 1:n
        ID = data(i, 1);
        ts_time = data(i, 2);
        te_time = data(i, 3);

        if ts_time < Ts1_time && te_time <= Te1_time && te_time > Ts1_time
            y1 = Ts1_time - ts_time;
            y2 = te_time - Ts1_time;
            y3 = 0;
            tag = 1;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts1_time && te_time <= Te1_time
            y1 = 0;
            y2 = te_time - ts_time;
            y3 = 0;
            tag = 2;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time < Ts1_time && te_time > Te1_time && te_time <= Te2_time
            y1 = Ts1_time - ts_time;
            y2 = Te1_time - Ts1_time;
            y3 = te_time - Te1_time;
            tag = 3;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts1_time && ts_time < Te1_time && te_time > Te1_time && te_time <= Te2_time
            y1 = 0;
            y2 = Te1_time - ts_time;
            y3 = te_time - Te1_time;
            tag = 4;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts1_time && ts_time < Te1_time && te_time > Te2_time
            y1 = 0;
            y2 = Te1_time - ts_time;
            y3 = Te2_time - Ts2_time;
            tag = 5;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time < Ts1_time && te_time > Te2_time
            y1 = Ts1_time - ts_time;
            y2 = Te1_time - Ts1_time;
            y3 = Te2_time - Ts2_time;
            tag = 6;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts2_time && te_time <= Te2_time
            y1 = 0;
            y2 = 0;
            y3 = te_time - ts_time;
            tag = 7;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts2_time && ts_time < Te2_time && te_time > Te2_time
            y1 = 0;
            y2 = 0;
            y3 = Te2_time - ts_time;
            tag = 8;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
    end
    
    % remove instances not related to W1 and W2
    IDs_delete = [];
    for i=1:n
        if new_data(i,5) == 0
            IDs_delete = [IDs_delete, i];
        end
    end
    new_data(IDs_delete,:) = [];      
end