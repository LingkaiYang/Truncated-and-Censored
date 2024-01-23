function new_data = data_extraction(data, W1, W2) 
    % data: original dataset, ID, start time and end time
    % W1: the start and end of the 1st window
    % W2: the start and end of the 2nd window
    Ts1 = W1(1);
    Te1 = W1(2);
    Ts2 = W2(1);
    Te2 = W2(2);
    
    % Ts1: the start of the 1st observation window, for example, "12/16/2012 18:33"
    % Te1: the end of the 1st observation window
    % Ts2: the start of the 2nd observation window
    % Te2: the end of the 2nd observation window
    
    format = 'MM/dd/yyyy HH:mm';
    Ts1_time = datetime(Ts1,'InputFormat',format);
    Te1_time = datetime(Te1,'InputFormat',format);
    Ts2_time = datetime(Ts2,'InputFormat',format);
    Te2_time = datetime(Te2,'InputFormat',format);
    
    [n,~] = size(data);
    new_data = zeros(n,5);

    for i = 1:n
        ID = str2double(data(i, 1));
        ts = data(i, 2);
        te = data(i, 3);
        ts_time = datetime(ts, 'InputFormat', format);
        te_time = datetime(te, 'InputFormat', format);

        if ts_time < Ts1_time && te_time <= Te1_time && te_time > Ts1_time
            y1 = hours(Ts1_time - ts_time);
            y2 = hours(te_time - Ts1_time);
            y3 = 0;
            tag = 1;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts1_time && te_time <= Te1_time
            y1 = 0;
            y2 = hours(te_time - ts_time);
            y3 = 0;
            tag = 2;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time < Ts1_time && te_time > Te1_time && te_time <= Te2_time
            y1 = hours(Ts1_time - ts_time);
            y2 = hours(Te1_time - Ts1_time);
            y3 = hours(te_time - Te1_time);
            tag = 3;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts1_time && ts_time < Te1_time && te_time > Te1_time && te_time <= Te2_time
            y1 = 0;
            y2 = hours(Te1_time - ts_time);
            y3 = hours(te_time - Te1_time);
            tag = 4;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts1_time && ts_time < Te1_time && te_time > Te2_time
            y1 = 0;
            y2 = hours(Te1_time - ts_time);
            y3 = hours(Te2_time - Ts2_time);
            tag = 5;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time < Ts1_time && te_time > Te2_time
            y1 = hours(Ts1_time - ts_time);
            y2 = hours(Te1_time - Ts1_time);
            y3 = hours(Te2_time - Ts2_time);
            tag = 6;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts2_time && te_time <= Te2_time
            y1 = 0;
            y2 = 0;
            y3 = hours(te_time - ts_time);
            tag = 7;
            new_data(i,:) = [ID, y1, y2, y3, tag];
        end
        if ts_time >= Ts2_time && ts_time < Te2_time && te_time > Te2_time
            y1 = 0;
            y2 = 0;
            y3 = hours(Te2_time - ts_time);
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