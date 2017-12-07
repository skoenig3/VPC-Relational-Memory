% analysis for alternative VPC task with no delay and long delays
% written by Seth Konig 7/2/15. Also for the RM project
%
clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';

%---Pre-Lesion---%
%combine all monkeys since don't have enough sets to run on individual
%monkeys
% pre_files = {'PW150710.2','PW150713.2','PW150714.2',...
%     'RR150710.2','RR150713.2','RR150714.2',...
%     'TO150710.2','TO150713.2','TO150714.2',...
%     'MF170414.2','MF170418.2','MF170419.2'};
% 
% sets = [7 9 11 ...
%     7 9 11 ...
%     8 10 12 ...
%     8 10 12];

pre_files = {'PW150710.2','PW150713.2','PW150714.2',...
    'RR150710.2','RR150713.2','RR150714.2',...
    'TO150710.2','TO150713.2','TO150714.2'};

%---Post-Lesion---%
post_files = {'TO170815.2','TO170817.2','TO170821.2',...
    'PW160531.2','PW160601.2',...
    'RR160622.2','RR160623.2','RR160624.2'};


samprate = 5/1000; % in secs
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
buffer = 100/samprate/1000; %100 ms buffer for filtering

%---Pre-Lesion Data---%
%AB
pre_prop_test = NaN(2,length(pre_files)); %row 1 short, row 2 long
pre_prop_test11 = NaN(2,length(pre_files)); % for first 250-1000 ms
pre_prop_test12 = NaN(2,length(pre_files)); % for first 1s - 5s
pre_test_transitions = NaN(2,length(pre_files));
pre_prop_test_time_short = zeros(length(pre_files),5000);
pre_prop_test_time_long  = zeros(length(pre_files),5000);
pre_valid_time_short = zeros(length(pre_files),5000);
pre_valid_time_long = zeros(length(pre_files),5000);

%AC
pre_prop_test2 = NaN(2,length(pre_files)); %row 1 short, row 2 long
pre_prop_test21 = NaN(2,length(pre_files)); % for first 250-1000 ms
pre_prop_test22 = NaN(2,length(pre_files)); % for first 1s - 5s
pre_test_transitions2 = NaN(2,length(pre_files));
pre_prop_test_time_short2 = zeros(length(pre_files),5000);
pre_prop_test_time_long2  = zeros(length(pre_files),5000);
pre_valid_time_short2 = zeros(length(pre_files),5000);
pre_valid_time_long2 = zeros(length(pre_files),5000);

for file = 1:length(pre_files)
    
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(10) '-fixation.mat']);
    
    
    pr_pl = NaN(2,60);
    trans= NaN(2,60);
    pr_pl11 = NaN(2,60);
    pr_pl12 = NaN(2,60);
    pr_pl21 = NaN(2,60);
    pr_pl22 = NaN(2,60);
    
    for trial = 2:3:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        %---Smooth Eye Data---%
        parsed_eyedat = preparse(eyedat{trial});
        smoothed_xy = [];
        for p = 1:length(parsed_eyedat)
            if any(~isnan(parsed_eyedat{p}(1,:)))
                
                %raw eye data
                x = parsed_eyedat{p}(1,:);
                y = parsed_eyedat{p}(2,:);
                
                %buffer before smoothing
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                
                %upsample
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                
                %low pass filtered velocity 30 Hz cutoff
                xss = filtfilt(flt,1,x);
                yss = filtfilt(flt,1,y);
                xss = xss(101:end-101); %remove buffer after filtering
                yss = yss(101:end-101); %remove buffer after filtering
                
                smoothed_xy =[smoothed_xy [xss; yss]];
            else
                smoothed_xy =[smoothed_xy parsed_eyedat{p}];
            end
        end
        x = smoothed_xy(1,:);
        y = smoothed_xy(2,:);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        
        if length(x) > 5000
            x = x(1:5000);
            y = y(1:5000);
        end
                        
        %only look at first 250-1000 ms or last 4 seconds
        x1 = x(250:1000);
        x2 = x(1001:5000);
        
        
    if img_pos(2,trial) == -12 %test image is on the left
            prl = x < -1;
            prl1 = sum(x1 < -1);
            prl2 = sum(x2 < -1);

        else %image is on the right
            prl = x > 1;
            prl1 = sum(x1 > 1);
            prl2 = sum(x2 > 1);
        end
        
        %number of transitions
        timesin = findgaps(find(prl));
        trans_count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                trans_count = trans_count+1;
            end
        end
        
        valid_time = x < -1 | x > 1;
        valid_time1 = (sum(x1 < -1) + sum(x1 > 1));
        valid_time2 = (sum(x2 < -1) + sum(x2 > 1));
        
        if trial <= 45 && trial >= 18 %short gaps ~10 secs
            pr_pl(1,trial) = sum(prl)/sum(valid_time);
            pr_pl11(1,trial) = prl1/valid_time1;
            pr_pl12(1,trial) = prl2/valid_time2;
            trans(1,trial) = trans_count;
            
            pre_prop_test_time_short(file,prl) = pre_prop_test_time_short(file,prl)+1;
            pre_valid_time_short(file,valid_time) = pre_valid_time_short(file,valid_time)+1;
            
        else %long gap > 2 minutes
            pr_pl(2,trial) = sum(prl)/sum(valid_time);
            pr_pl11(2,trial) = prl1/valid_time1;
            pr_pl12(2,trial) = prl2/valid_time2;
            trans(2,trial) = trans_count;
            
            pre_prop_test_time_long(file,prl) = pre_prop_test_time_long (file,prl)+1;
            pre_valid_time_long(file,valid_time) = pre_valid_time_long(file,valid_time)+1;
        end
        
    end
    
    for trial = 3:3:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        %---Smooth Eye Data---%
        parsed_eyedat = preparse(eyedat{trial});
        smoothed_xy = [];
        for p = 1:length(parsed_eyedat)
            if any(~isnan(parsed_eyedat{p}(1,:)))
                
                %raw eye data
                x = parsed_eyedat{p}(1,:);
                y = parsed_eyedat{p}(2,:);
                
                %buffer before smoothing
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                
                %upsample
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                
                %low pass filtered velocity 30 Hz cutoff
                xss = filtfilt(flt,1,x);
                yss = filtfilt(flt,1,y);
                xss = xss(101:end-101); %remove buffer after filtering
                yss = yss(101:end-101); %remove buffer after filtering
                
                smoothed_xy =[smoothed_xy [xss; yss]];
            else
                smoothed_xy =[smoothed_xy parsed_eyedat{p}];
            end
        end
        x = smoothed_xy(1,:);
        y = smoothed_xy(2,:);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        
        if length(x) > 5000
            x = x(1:5000);
            y = y(1:5000);
        end
                
        %only look at first 250-1000 ms or last 4 seconds
        x1 = x(250:1000);
        x2 = x(1001:5000);
        
        if img_pos(2,trial) == -12 %test image is on the left
            prl = x < -1;
            prl1 = sum(x1 < -1);
            prl2 = sum(x2 < -1);
            
        else %image is on the right
            prl = x > 1;
            prl1 = sum(x1 > 1);
            prl2 = sum(x2 > 1);
        end
        
        %number of transitions
        timesin = findgaps(find(prl));
        trans_count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                trans_count = trans_count+1;
            end
        end
        
        valid_time = x < -1 | x > 1;
        valid_time1 = (sum(x1 < -1) + sum(x1 > 1));
        valid_time2 = (sum(x2 < -1) + sum(x2 > 1));
        
        if trial <= 45 && trial >= 18 %short gaps ~10 secs
            pr_pl(1,trial) = sum(prl)/sum(valid_time);
            pr_pl21(1,trial) = prl1/valid_time1;
            pr_pl22(1,trial) = prl2/valid_time2;
            trans(1,trial) = trans_count;
            
            pre_prop_test_time_short2(file,prl) = pre_prop_test_time_short(file,prl)+1;
            pre_valid_time_short2(file,valid_time) = pre_valid_time_short(file,valid_time)+1;
            
        else %long gap > 2 minutes
            pr_pl(2,trial) = sum(prl)/sum(valid_time);
            pr_pl21(2,trial) = prl1/valid_time1;
            pr_pl22(2,trial) = prl2/valid_time2;
            trans(2,trial) = trans_count;
            
            pre_prop_test_time_long2(file,prl) = pre_prop_test_time_long (file,prl)+1;
            pre_valid_time_long2(file,valid_time) = pre_valid_time_long(file,valid_time)+1;
        end
        
    end
    
    %reformat for whole trial
    prplAB = pr_pl(:,2:3:end);
    prplAC = pr_pl(:,3:3:end);
    transAB = trans(:,2:3:end);
    transAC = trans(:,3:3:end);


    %AB
    pre_prop_test(1,file) = nanmean(prplAB(1,:));
    pre_prop_test(2,file) = nanmean(prplAB(2,:));
    pre_test_transitions(1,file) = nanmean(transAB(1,:));
    pre_test_transitions(2,file) = nanmean(transAB(2,:));
    
    %AC
    pre_prop_test2(1,file) = nanmean(prplAC(1,:));
    pre_prop_test2(2,file) = nanmean(prplAC(2,:));
    pre_test_transitions2(1,file) = nanmean(transAC(1,:));
    pre_test_transitions2(2,file) = nanmean(transAC(2,:));
    
    
    %for various time window
    pre_prop_test11(1,file) = nanmean(pr_pl11(1,:));
    pre_prop_test11(2,file) = nanmean(pr_pl11(2,:));
    pre_prop_test12(1,file) = nanmean(pr_pl12(1,:));
    pre_prop_test12(2,file) = nanmean(pr_pl12(2,:));
    pre_prop_test21(1,file) = nanmean(pr_pl21(1,:));
    pre_prop_test21(2,file) = nanmean(pr_pl21(2,:));
    pre_prop_test22(1,file) = nanmean(pr_pl22(1,:));
    pre_prop_test22(2,file) = nanmean(pr_pl22(2,:));
    
end

%% Recalculate proportion of time on novel stimulus by dividing by total time looking at either stimulus
pre_prop_time_short = zeros(size(pre_prop_test_time_short));
for file = 1:size(pre_prop_test_time_short,1)
    valid_time = pre_valid_time_short(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_short(file,:) = pre_prop_test_time_short(file,:)./valid_time;
end

pre_prop_time_long = zeros(size(pre_prop_test_time_long));
for file = 1:size(pre_prop_test_time_long,1)
    valid_time = pre_valid_time_long(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_long(file,:) = pre_prop_test_time_long(file,:)./valid_time;
end

pre_prop_time_short2 = zeros(size(pre_prop_test_time_short2));
for file = 1:size(pre_prop_test_time_short2,1)
    valid_time = pre_valid_time_short2(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_short2(file,:) = pre_prop_test_time_short2(file,:)./valid_time;
end

pre_prop_time_long2 = zeros(size(pre_prop_test_time_long2));
for file = 1:size(pre_prop_test_time_long2,1)
    valid_time = pre_valid_time_long2(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_long2(file,:) = pre_prop_test_time_long2(file,:)./valid_time;
end
%%
%---Pre-Lesion Data---%
%AB
post_prop_test = NaN(2,length(post_files)); %row 1 short, row 2 long
post_prop_test11 = NaN(2,length(post_files)); % for first 250-1000 ms
post_prop_test12 = NaN(2,length(post_files)); % for first 1s - 5s
post_test_transitions = NaN(2,length(post_files));
post_prop_test_time_short = zeros(length(post_files),5000);
post_prop_test_time_long  = zeros(length(post_files),5000);
post_valid_time_short = zeros(length(post_files),5000);
post_valid_time_long = zeros(length(post_files),5000);

%AC
post_prop_test2 = NaN(2,length(post_files)); %row 1 short, row 2 long
post_prop_test21 = NaN(2,length(post_files)); % for first 250-1000 ms
post_prop_test22 = NaN(2,length(post_files)); % for first 1s - 5s
post_test_transitions2 = NaN(2,length(post_files));
post_prop_test_time_short2 = zeros(length(post_files),5000);
post_prop_test_time_long2  = zeros(length(post_files),5000);
post_valid_time_short2 = zeros(length(post_files),5000);
post_valid_time_long2 = zeros(length(post_files),5000);

for file = 1:length(post_files)
    
    load([data_dir post_files{file}(1:8) '_' post_files{file}(10) '-fixation.mat']);
    
    
    pr_pl = NaN(2,60);
    trans= NaN(2,60);
    pr_pl11 = NaN(2,60);
    pr_pl12 = NaN(2,60);
    pr_pl21 = NaN(2,60);
    pr_pl22 = NaN(2,60);
    
    for trial = 2:3:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        %---Smooth Eye Data---%
        parsed_eyedat = preparse(eyedat{trial});
        smoothed_xy = [];
        for p = 1:length(parsed_eyedat)
            if any(~isnan(parsed_eyedat{p}(1,:)))
                
                %raw eye data
                x = parsed_eyedat{p}(1,:);
                y = parsed_eyedat{p}(2,:);
                
                %buffer before smoothing
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                
                %upsample
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                
                %low pass filtered velocity 30 Hz cutoff
                xss = filtfilt(flt,1,x);
                yss = filtfilt(flt,1,y);
                xss = xss(101:end-101); %remove buffer after filtering
                yss = yss(101:end-101); %remove buffer after filtering
                
                smoothed_xy =[smoothed_xy [xss; yss]];
            else
                smoothed_xy =[smoothed_xy parsed_eyedat{p}];
            end
        end
        x = smoothed_xy(1,:);
        y = smoothed_xy(2,:);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        
        if length(x) > 5000
            x = x(1:5000);
            y = y(1:5000);
        end
                        
        %only look at first 250-1000 ms or last 4 seconds
        x1 = x(250:1000);
        x2 = x(1001:5000);
        
        
    if img_pos(2,trial) == -12 %test image is on the left
            prl = x < -1;
            prl1 = sum(x1 < -1);
            prl2 = sum(x2 < -1);

        else %image is on the right
            prl = x > 1;
            prl1 = sum(x1 > 1);
            prl2 = sum(x2 > 1);
        end
        
        %number of transitions
        timesin = findgaps(find(prl));
        trans_count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                trans_count = trans_count+1;
            end
        end
        
        valid_time = x < -1 | x > 1;
        valid_time1 = (sum(x1 < -1) + sum(x1 > 1));
        valid_time2 = (sum(x2 < -1) + sum(x2 > 1));
        
        if trial <= 45 && trial >= 18 %short gaps ~10 secs
            pr_pl(1,trial) = sum(prl)/sum(valid_time);
            pr_pl11(1,trial) = prl1/valid_time1;
            pr_pl12(1,trial) = prl2/valid_time2;
            trans(1,trial) = trans_count;
            
            post_prop_test_time_short(file,prl) = post_prop_test_time_short(file,prl)+1;
            post_valid_time_short(file,valid_time) = post_valid_time_short(file,valid_time)+1;
            
        else %long gap > 2 minutes
            pr_pl(2,trial) = sum(prl)/sum(valid_time);
            pr_pl11(2,trial) = prl1/valid_time1;
            pr_pl12(2,trial) = prl2/valid_time2;
            trans(2,trial) = trans_count;
            
            post_prop_test_time_long(file,prl) = post_prop_test_time_long (file,prl)+1;
            post_valid_time_long(file,valid_time) = post_valid_time_long(file,valid_time)+1;
        end
        
    end
    
    for trial = 3:3:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        %---Smooth Eye Data---%
        parsed_eyedat = preparse(eyedat{trial});
        smoothed_xy = [];
        for p = 1:length(parsed_eyedat)
            if any(~isnan(parsed_eyedat{p}(1,:)))
                
                %raw eye data
                x = parsed_eyedat{p}(1,:);
                y = parsed_eyedat{p}(2,:);
                
                %buffer before smoothing
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                
                %upsample
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                
                %low pass filtered velocity 30 Hz cutoff
                xss = filtfilt(flt,1,x);
                yss = filtfilt(flt,1,y);
                xss = xss(101:end-101); %remove buffer after filtering
                yss = yss(101:end-101); %remove buffer after filtering
                
                smoothed_xy =[smoothed_xy [xss; yss]];
            else
                smoothed_xy =[smoothed_xy parsed_eyedat{p}];
            end
        end
        x = smoothed_xy(1,:);
        y = smoothed_xy(2,:);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        
        if length(x) > 5000
            x = x(1:5000);
            y = y(1:5000);
        end
                
        %only look at first 250-1000 ms or last 4 seconds
        x1 = x(250:1000);
        x2 = x(1001:5000);
        
        if img_pos(2,trial) == -12 %test image is on the left
            prl = x < -1;
            prl1 = sum(x1 < -1);
            prl2 = sum(x2 < -1);
            
        else %image is on the right
            prl = x > 1;
            prl1 = sum(x1 > 1);
            prl2 = sum(x2 > 1);
        end
        
        %number of transitions
        timesin = findgaps(find(prl));
        trans_count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                trans_count = trans_count+1;
            end
        end
        
        valid_time = x < -1 | x > 1;
        valid_time1 = (sum(x1 < -1) + sum(x1 > 1));
        valid_time2 = (sum(x2 < -1) + sum(x2 > 1));
        
        if trial <= 45 && trial >= 18 %short gaps ~10 secs
            pr_pl(1,trial) = sum(prl)/sum(valid_time);
            pr_pl21(1,trial) = prl1/valid_time1;
            pr_pl22(1,trial) = prl2/valid_time2;
            trans(1,trial) = trans_count;
            
            post_prop_test_time_short2(file,prl) = post_prop_test_time_short(file,prl)+1;
            post_valid_time_short2(file,valid_time) = post_valid_time_short(file,valid_time)+1;
            
        else %long gap > 2 minutes
            pr_pl(2,trial) = sum(prl)/sum(valid_time);
            pr_pl21(2,trial) = prl1/valid_time1;
            pr_pl22(2,trial) = prl2/valid_time2;
            trans(2,trial) = trans_count;
            
            post_prop_test_time_long2(file,prl) = post_prop_test_time_long (file,prl)+1;
            post_valid_time_long2(file,valid_time) = post_valid_time_long(file,valid_time)+1;
        end
        
    end
    
    %reformat for whole trial
    prplAB = pr_pl(:,2:3:end);
    prplAC = pr_pl(:,3:3:end);
    transAB = trans(:,2:3:end);
    transAC = trans(:,3:3:end);


    %AB
    post_prop_test(1,file) = nanmean(prplAB(1,:));
    post_prop_test(2,file) = nanmean(prplAB(2,:));
    post_test_transitions(1,file) = nanmean(transAB(1,:));
    post_test_transitions(2,file) = nanmean(transAB(2,:));
    
    %AC
    post_prop_test2(1,file) = nanmean(prplAC(1,:));
    post_prop_test2(2,file) = nanmean(prplAC(2,:));
    post_test_transitions2(1,file) = nanmean(transAC(1,:));
    post_test_transitions2(2,file) = nanmean(transAC(2,:));
    
    
    %for various time window
    post_prop_test11(1,file) = nanmean(pr_pl11(1,:));
    post_prop_test11(2,file) = nanmean(pr_pl11(2,:));
    post_prop_test12(1,file) = nanmean(pr_pl12(1,:));
    post_prop_test12(2,file) = nanmean(pr_pl12(2,:));
    post_prop_test21(1,file) = nanmean(pr_pl21(1,:));
    post_prop_test21(2,file) = nanmean(pr_pl21(2,:));
    post_prop_test22(1,file) = nanmean(pr_pl22(1,:));
    post_prop_test22(2,file) = nanmean(pr_pl22(2,:));
    
end

%% Recalculate proportion of time on novel stimulus by dividing by total time looking at either stimulus
post_prop_time_short = zeros(size(post_prop_test_time_short));
for file = 1:size(post_prop_test_time_short,1)
    valid_time = post_valid_time_short(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_short(file,:) = post_prop_test_time_short(file,:)./valid_time;
end

post_prop_time_long = zeros(size(post_prop_test_time_long));
for file = 1:size(post_prop_test_time_long,1)
    valid_time = post_valid_time_long(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_long(file,:) = post_prop_test_time_long(file,:)./valid_time;
end

post_prop_time_short2 = zeros(size(post_prop_test_time_short2));
for file = 1:size(post_prop_test_time_short2,1)
    valid_time = post_valid_time_short2(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_short2(file,:) = post_prop_test_time_short2(file,:)./valid_time;
end

post_prop_time_long2 = zeros(size(post_prop_test_time_long2));
for file = 1:size(post_prop_test_time_long2,1)
    valid_time = post_valid_time_long2(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_long2(file,:) = post_prop_test_time_long2(file,:)./valid_time;
end

%%


figure
subplot(1,2,1)
hold on
dofill(1:5000,100*pre_prop_time_short/1000,'blue',1,50)
dofill(1:5000,100*post_prop_time_short/1000,'red',1,50)
dofill(1:5000,100*pre_prop_time_long/1000,'green',1,50)
dofill(1:5000,100*post_prop_time_long/1000,'magenta',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
title(['AB'])
box off
legend({'Pre-Short','Post-Short','Pre-Long','Post-Long'})

subplot(1,2,2)
hold on
dofill(1:5000,100*pre_prop_time_short2/1000,'blue',1,50)
dofill(1:5000,100*pre_prop_time_long2/1000,'red',1,50)
dofill(1:5000,100*post_prop_time_short2/1000,'green',1,50)
dofill(1:5000,100*post_prop_time_long2/1000,'magenta',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
title(['AC'])
box off
legend({'Pre-Short','Post-Short','Pre-Long','Post-Long'})

subtitle('Group VPC2SCW')


%%


figure
subplot(2,2,1)
errorb(100*[mean(pre_prop_test,2) mean(pre_prop_test2,2)],...
    100*[std(pre_prop_test,[],2)./sqrt(length(pre_files))...
    std(pre_prop_test2,[],2)./sqrt(length(pre_files))])
hold on
plot([0.5 2.5],[50 50],'k--')
hold off
box off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
ylabel('% of time spent looking at Novel Stimulus')
yl = ylim;
ylim([40 100])
legend('AB','AC')
title('Whole Trial')


subplot(2,2,2)
errorb([mean(pre_test_transitions,2) mean(pre_test_transitions2,2)],...
    [std(pre_test_transitions,[],2)./sqrt(length(pre_files))...
    std(pre_test_transitions2,[],2)./sqrt(length(pre_files))])
box off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
ylabel('# of Transitions in and out of Novel Stimulus')
title('Number of Transitions')


subplot(2,2,3)
errorb(100*[mean(pre_prop_test11,2) mean(pre_prop_test21,2)],...
    100*[std(pre_prop_test11,[],2)./sqrt(length(pre_files))...
    std(pre_prop_test21,[],2)./sqrt(length(pre_files))])
box off
hold on
plot([0.5 2.5],[50 50],'k--')
hold off
ylim([40 100])
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
ylabel('% of time spent looking at Novel Stimulus')
title('250-1000 ms')


subplot(2,2,4)
errorb(100*[mean(pre_prop_test12,2) mean(pre_prop_test22,2)],...
    100*[std(pre_prop_test22,[],2)./sqrt(length(pre_files))...
    std(pre_prop_test22,[],2)./sqrt(length(pre_files))])
hold on
plot([0.5 2.5],[50 50],'k--')
hold off
box off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
yl = ylim;
ylim([40 100])
legend('Pre','Post')
title('1000-5000 ms')
subtitle(['VPC Black/White Wide Apart: Novelty Preference Over Time: All Monkeys Pre-Lesion'])

