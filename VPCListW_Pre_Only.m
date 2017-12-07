% analysis for alternative VPC task with no delay and long delays
% written by Seth Konig 7/2/15. Also for the RM project

clar

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';

%all_monkeys
pre_files = {'PW150715.2','PW150716.2','PW150717.2',...
    'RR150715.2','RR150716.2','RR150720.2',...
    'TO150715.2','TO150716.2','TO150717.2'};

post_files = {'PW160602.2','PW160603.2','PW160606.2',...
    'RR160627.2','RR160629.2','RR160630.2',...
    'TO170822.2','TO170823.2','TO170828.2'};

samprate = 5/1000; % in secs
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
buffer = 100/samprate/1000; %100 ms buffer for filtering
%%
%---Pre-Lesion Data---%
pre_prop_test = NaN(length(pre_files),20); %row 1 short, row 2 long
pre_prop_test_time_short = zeros(length(pre_files),5000); %<= 12
pre_prop_test_time_medium = zeros(length(pre_files),5000); %remaining
pre_prop_test_time_long  = zeros(length(pre_files),5000); %> 24
pre_valid_time_short = zeros(length(pre_files),5000);
pre_valid_time_medium = zeros(length(pre_files),5000);
pre_valid_time_long = zeros(length(pre_files),5000);

for file = 1:length(pre_files)
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(10) '-fixation.mat']);
    
    pl_pr = NaN(1,20);
    del = NaN(1,20);
    
    for imnum = 1:20
        [row,img_trial] = find(imgnum == imnum);
        test_row = row(end);
        test_img_trial = img_trial(end);
        
        x = eyedat{test_img_trial}(1,:);
        y = eyedat{test_img_trial}(2,:);
        
        %---Smooth Eye Data---%
        parsed_eyedat = preparse(eyedat{test_img_trial});
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
        
        if test_row == 2 %then novel image on left
            prl = x < -1;
        else %image is on the right
            prl = x > 1;
        end
        valid_time = x < -1 | x > 1;
        
        pl_pr(imnum) = sum(prl(250:1000))/sum(valid_time(250:1000));
        del(imnum) = img_trial(end)-img_trial(1)-1; %number of intervening trials
        
        if del(imnum) <= 12
            pre_prop_test_time_short(file,prl) = pre_prop_test_time_short(file,prl)+1;
            pre_valid_time_short(file,valid_time) = pre_valid_time_short(file,valid_time)+1;
        elseif del(imnum) > 24
            pre_prop_test_time_long(file,prl) = pre_prop_test_time_long(file,prl)+1;
            pre_valid_time_long(file,valid_time) = pre_valid_time_long(file,valid_time)+1;
        else
            pre_prop_test_time_medium(file,prl) = pre_prop_test_time_medium(file,prl)+1;
            pre_valid_time_medium(file,valid_time) = pre_valid_time_medium(file,valid_time)+1;
        end
    end
    
    [~,sort_ind] = sort(del); %delays are not perfectly in order for some reason but seem ok overall
    
    pre_prop_test(file,:) = pl_pr(sort_ind);
end
%% Recalculate proportion of time on novel stimulus by dividing by total time looking at either stimulus

pre_prop_time_short = zeros(size(pre_prop_test_time_short));
for file = 1:size(pre_prop_test_time_short,1)
    valid_time = pre_valid_time_short(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_short(file,:) = pre_prop_test_time_short(file,:)./valid_time;
end

pre_prop_time_medium = zeros(size(pre_prop_test_time_medium));
for file = 1:size(pre_prop_test_time_medium,1)
    valid_time = pre_valid_time_medium(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_medium(file,:) = pre_prop_test_time_medium(file,:)./valid_time;
end

pre_prop_time_long = zeros(size(pre_prop_test_time_long));
for file = 1:size(pre_prop_test_time_long,1)
    valid_time = pre_valid_time_long(file,:);
    valid_time(valid_time < 5) = 5;
    pre_prop_time_long(file,:) = pre_prop_test_time_long(file,:)./valid_time;
end
%%
%---Pre-Lesion Data---%
post_prop_test = NaN(length(post_files),20); %row 1 short, row 2 long
post_prop_test_time_short = zeros(length(post_files),5000); %<= 12
post_prop_test_time_medium = zeros(length(post_files),5000); %remaining
post_prop_test_time_long  = zeros(length(post_files),5000); %> 24
post_valid_time_short = zeros(length(post_files),5000);
post_valid_time_medium = zeros(length(post_files),5000);
post_valid_time_long = zeros(length(post_files),5000);

for file = 1:length(post_files)
    load([data_dir post_files{file}(1:8) '_' post_files{file}(10) '-fixation.mat']);
    
    pl_pr = NaN(1,20);
    del = NaN(1,20);
    
    for imnum = 1:20
        [row,img_trial] = find(imgnum == imnum);
        test_row = row(end);
        test_img_trial = img_trial(end);
        
        x = eyedat{test_img_trial}(1,:);
        y = eyedat{test_img_trial}(2,:);
        
        %---Smooth Eye Data---%
        parsed_eyedat = preparse(eyedat{test_img_trial});
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
        
        if test_row == 2 %then novel image on left
            prl = x < -1;
        else %image is on the right
            prl = x > 1;
        end
        valid_time = x < -1 | x > 1;
        
        pl_pr(imnum) = sum(prl(250:1000))/sum(valid_time(250:1000));
        del(imnum) = img_trial(end)-img_trial(1)-1; %number of intervening trials
        
        if del(imnum) <= 12
            post_prop_test_time_short(file,prl) = post_prop_test_time_short(file,prl)+1;
            post_valid_time_short(file,valid_time) = post_valid_time_short(file,valid_time)+1;
        elseif del(imnum) > 24
            post_prop_test_time_long(file,prl) = post_prop_test_time_long(file,prl)+1;
            post_valid_time_long(file,valid_time) = post_valid_time_long(file,valid_time)+1;
        else
            post_prop_test_time_medium(file,prl) = post_prop_test_time_medium(file,prl)+1;
            post_valid_time_medium(file,valid_time) = post_valid_time_medium(file,valid_time)+1;
        end
    end
    
    [~,sort_ind] = sort(del); %delays are not perfectly in order for some reason but seem ok overall
    
    post_prop_test(file,:) = pl_pr(sort_ind);
end
%% Recalculate proportion of time on novel stimulus by dividing by total time looking at either stimulus

post_prop_time_short = zeros(size(post_prop_test_time_short));
for file = 1:size(post_prop_test_time_short,1)
    valid_time = post_valid_time_short(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_short(file,:) = post_prop_test_time_short(file,:)./valid_time;
end

post_prop_time_medium = zeros(size(post_prop_test_time_medium));
for file = 1:size(post_prop_test_time_medium,1)
    valid_time = post_valid_time_medium(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_medium(file,:) = post_prop_test_time_medium(file,:)./valid_time;
end

post_prop_time_long = zeros(size(post_prop_test_time_long));
for file = 1:size(post_prop_test_time_long,1)
    valid_time = post_valid_time_long(file,:);
    valid_time(valid_time < 5) = 5;
    post_prop_time_long(file,:) = post_prop_test_time_long(file,:)./valid_time;
end

%%
figure
subplot(1,2,1)
hold on
dofill(1:5000,100*pre_prop_time_short/1000,'blue',1,50)
dofill(1:5000,100*pre_prop_time_medium/1000,'green',1,50)
dofill(1:5000,100*pre_prop_time_long/1000,'red',1,50)
dofill(1:5000,100*post_prop_time_short/1000,'magenta',1,50)
dofill(1:5000,100*post_prop_time_medium/1000,'black',1,50)
dofill(1:5000,100*post_prop_time_long/1000,'yellow',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
title(['Novelty Preference Over Time'])
box off
legend({'Pre-Short','Pre-Medium','Pre-Long','Post-Short','Post-Medium','Post-Long'})

subplot(1,2,2)
hold on
errorb(1:20,100*mean(pre_prop_test),100*std(pre_prop_test)./sqrt(length(pre_files)))
errorb(1:20+0.5,100*mean(post_prop_test),100*std(post_prop_test)./sqrt(length(post_files)),'color','r')
plot([0 21],[50 50],'k--')
hold off
xlim([0 21])
xlabel('# of Interleaving Stimuli')
ylabel('Novelty Preference')
ylim([30 100])
title('Average Novelty 250-1000 ms')
box off
legend('Pre','Post')

subtitle('All Monkeys List VPC')
