%written Seth Konig 6/9/16
% requires that files be processed in Batch_VPCRM

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';


%---Vivian---%
pre_files = {'PW150609.2','PW150610.2','PW150611.2','PW150612.2',...
    'PW150615.2','PW150616.2','PW150617.2','PW150618.2',...
    'PW150619.2','PW150622.2','PW150623.2','PW150624.2',...
    'PW150625.2','PW150626.2','PW150629.2','PW150630.2'};

post_files = {'PW160425.2','PW160426.2','PW160427.2','PW160428.2','PW160429.2',...
    'PW160502.2','PW160503.2','PW160504.2','PW160505.2','PW160506.2',...
    'PW160509.2','PW160513.2','PW160516.2','PW160518.2','PW160519.2'};

%---Red---%
% pre_files  = {'RR150609.2','RR150610.2','RR150611.2','RR150612.2',...
%              'RR150615.2','RR150616.2','RR150617.2','RR150618.2',...
%              'RR150619.2','RR150622.2','RR150623.2','RR150624.2',...
%              'RR150625.2','RR150626.2','RR150629.2','RR150630.2'};
% 
% post_files = {'RR160523.2','RR160526.2','RR160531.2','RR160601.2',...
%              'RR160602.2','RR160603.2','RR160606.2','RR160607.2',...
%              'RR160608.2','RR160609.2','RR160610.2','RR160613.2',...
%              'RR160614.2','RR160615.2'};

%---Get Data for Familiariztion Phase (Same Novel image on both sides)---%
pre_novel_transitions = NaN(2,length(pre_files));
pre_novel_propleft_propright = NaN(2,length(pre_files)); %for whole session
pre_novel_propleft_propright1 = NaN(2,length(pre_files));%for first ~1 second (250 ms -1000 ms)
for file = 1:length(pre_files)
    pl_pr = NaN(2,40);
    pl_pr1 = NaN(2,40);
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(10) '-fixation.mat']);
    
    trans1 = []; %left
    trans2 = []; %right
    for trial = 1:2:length(eyedat); %every other trial is novel/fam phase
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        img_on = round(per(trial).alltim(per(trial).allval == 23)/5);
     
        badx = find(x < -8.75 | x > 8.75); 
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
               
        pl_pr(1,trial) = sum(x < -1)/(sum(x < -1) + sum(x > 1));
        pl_pr(2,trial) = sum(x >  1)/(sum(x < -1) + sum(x > 1));
        
        %only look at first 250-1000 ms
        x1 = x(50:200);
        pl_pr1(1,trial) = sum(x1 < -1)/(sum(x1 < -1) + sum(x1 > 1));
        pl_pr1(2,trial) = sum(x1 >  1)/(sum(x1 < -1) + sum(x1 > 1));
        
        timesin = findgaps(find(x < -1));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans1 = [trans1 count];
        
        timesin = findgaps(find(x > 1));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans2 = [trans2 count];
        
    end
    pre_novel_propleft_propright(:,file) = nanmean(pl_pr,2);
    pre_novel_propleft_propright1(:,file) = nanmean(pl_pr1,2);
    pre_novel_transitions(1,file) = mean(trans1); %left
    pre_novel_transitions(2,file) = mean(trans2); %right
end

%bayesian inference to get posterior usually very minor
side_bias_correction = 2*pre_novel_propleft_propright./[sum(pre_novel_propleft_propright); sum(pre_novel_propleft_propright)];
side_bias_correction1 = ones(2,length(pre_files));%2*pre_novel_propleft_propright1./[sum(pre_novel_propleft_propright1); sum(pre_novel_propleft_propright1)];

%---Get Data for Test Phase (1 Novel 1 repeat)---%
pre_prop_test = NaN(2,length(pre_files)); %row 1 short, row 2 long
pre_prop_test1 = NaN(2,length(pre_files)); % for first 250-1000 ms
pre_test_transitions = NaN(2,length(pre_files));
pre_prop_test_time_short = zeros(length(pre_files),1000);
pre_prop_test_time_long  = zeros(length(pre_files),1000);

for file = 1:length(pre_files)
    pr_pl = NaN(2,40);
    pr_pl1 = NaN(2,40);
    load([data_dir pre_files{file}(1:8) '_' pre_files{file}(10) '-fixation.mat']);
    
    trans_short = [];
    trans_long = [];
    for trial = 2:2:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        badx = find(x < -8.75 | x > 8.75);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
        
         %only look at first 250-1000 ms
        x1 = x(50:200);
       
        
        if img_pos(2,trial) == -5 %test image is on the left
           
            timesin = findgaps(find(x < -1));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            
            if trial <= 30 && trial >= 12 %short gaps ~10 secs
                pr_pl(1,trial) = sum(x < -1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(1,trial) = sum(x1 < -1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                pre_prop_test_time_short(file,x < -1) = pre_prop_test_time_short(file,x < -1)+1;
                trans_short = [trans_short count];
            else %long gaps ~ 5 minutes
                pr_pl(2,trial) = sum(x < -1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(2,trial) = sum(x1 < -1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                pre_prop_test_time_long (file,x < -1) = pre_prop_test_time_long (file,x < -1)+1;
                trans_long = [trans_long count];
            end
            
        else %test image is on the right

            timesin = findgaps(find(x > 1));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            
            if trial <= 30 && trial >= 12
                pr_pl(1,trial) = sum(x > 1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(1,trial) = sum(x1 > 1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                pre_prop_test_time_short(file,x > 1) = pre_prop_test_time_short(file,x > 1)+1;
                trans_short = [trans_short count];
            else
                pr_pl(2,trial) = sum(x > 1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(2,trial) = sum(x1 > 1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                pre_prop_test_time_long (file,x > 1) = pre_prop_test_time_long (file,x > 1)+1;
                trans_long = [trans_long count];
            end
            
        end
    end
    
    pre_prop_test(1,file) = nanmean(pr_pl(1,:));
    pre_prop_test(2,file) = nanmean(pr_pl(2,:));
    
    pre_prop_test1(1,file) = nanmean(pr_pl1(1,:));
    pre_prop_test1(2,file) = nanmean(pr_pl1(2,:));
    
    pre_test_transitions(1,file) = mean(trans_short);
    pre_test_transitions(2,file) = mean(trans_long);
end

%---Get Data for Familiariztion Phase (Same Novel image on both sides)---%
post_novel_transitions = NaN(2,length(post_files));
post_novel_propleft_propright = NaN(2,length(post_files)); %for whole session
post_novel_propleft_propright1 = NaN(2,length(post_files));%for first ~1 second (250 ms -1000 ms)
for file = 1:length(post_files)
    pl_pr = NaN(2,40);
    pl_pr1 = NaN(2,40);
    load([data_dir post_files{file}(1:8) '_' post_files{file}(10) '-fixation.mat']);
    
    trans1 = []; %left
    trans2 = []; %right
    for trial = 1:2:length(eyedat); %every other trial is novel/fam phase
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        img_on = round(per(trial).alltim(per(trial).allval == 23)/5);
     
        badx = find(x < -8.75 | x > 8.75); 
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
               
        pl_pr(1,trial) = sum(x < -1)/(sum(x < -1) + sum(x > 1));
        pl_pr(2,trial) = sum(x >  1)/(sum(x < -1) + sum(x > 1));
        
        %only look at first 250-1000 ms
        x1 = x(50:200);
        pl_pr1(1,trial) = sum(x1 < -1)/(sum(x1 < -1) + sum(x1 > 1));
        pl_pr1(2,trial) = sum(x1 >  1)/(sum(x1 < -1) + sum(x1 > 1));
        
        timesin = findgaps(find(x < -1));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans1 = [trans1 count];
        
        timesin = findgaps(find(x > 1));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans2 = [trans2 count];
        
    end
    post_novel_propleft_propright(:,file) = nanmean(pl_pr,2);
    post_novel_propleft_propright1(:,file) = nanmean(pl_pr1,2);
    post_novel_transitions(1,file) = mean(trans1); %left
    post_novel_transitions(2,file) = mean(trans2); %right
end

%bayesian inference to get posterior usually very minor
side_bias_correction = 2*post_novel_propleft_propright./[sum(post_novel_propleft_propright); sum(post_novel_propleft_propright)];
side_bias_correction1 = ones(2,length(post_files));%2*pre_novel_propleft_propright1./[sum(pre_novel_propleft_propright1); sum(pre_novel_propleft_propright1)];

%---Get Data for Test Phase (1 Novel 1 repeat)---%
post_prop_test = NaN(2,length(post_files)); %row 1 short, row 2 long
post_prop_test1 = NaN(2,length(post_files)); % for first 250-1000 ms
post_test_transitions = NaN(2,length(post_files));
post_prop_test_time_short = zeros(length(post_files),1000);
post_prop_test_time_long  = zeros(length(post_files),1000);

for file = 1:length(post_files)
    pr_pl = NaN(2,40);
    pr_pl1 = NaN(2,40);
    load([data_dir post_files{file}(1:8) '_' post_files{file}(10) '-fixation.mat']);
    
    trans_short = [];
    trans_long = [];
    for trial = 2:2:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        badx = find(x < -8.75 | x > 8.75);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = NaN;
        y(badind) = NaN;
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
        
         %only look at first 250-1000 ms
        x1 = x(50:200);
       
        
        if img_pos(2,trial) == -5 %test image is on the left
           
            timesin = findgaps(find(x < -1));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            
            if trial <= 30 && trial >= 12 %short gaps ~10 secs
                pr_pl(1,trial) = sum(x < -1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(1,trial) = sum(x1 < -1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                post_prop_test_time_short(file,x < -1) = post_prop_test_time_short(file,x < -1)+1;
                trans_short = [trans_short count];
            else %long gaps ~ 5 minutes
                pr_pl(2,trial) = sum(x < -1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(2,trial) = sum(x1 < -1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                post_prop_test_time_long (file,x < -1) = post_prop_test_time_long (file,x < -1)+1;
                trans_long = [trans_long count];
            end
            
        else %test image is on the right

            timesin = findgaps(find(x > 1));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            
            if trial <= 30 && trial >= 12
                pr_pl(1,trial) = sum(x > 1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(1,trial) = sum(x1 > 1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                post_prop_test_time_short(file,x > 1) = post_prop_test_time_short(file,x > 1)+1;
                trans_short = [trans_short count];
            else
                pr_pl(2,trial) = sum(x > 1)/(sum(x < -1) + sum(x > 1))*side_bias_correction(1,file);
                pr_pl1(2,trial) = sum(x1 > 1)/(sum(x1 < -1) + sum(x1 > 1))*side_bias_correction1(1,file);
                post_prop_test_time_long (file,x > 1) = post_prop_test_time_long (file,x > 1)+1;
                trans_long = [trans_long count];
            end
            
        end
    end
    
    post_prop_test(1,file) = nanmean(pr_pl(1,:));
    post_prop_test(2,file) = nanmean(pr_pl(2,:));
    
    post_prop_test1(1,file) = nanmean(pr_pl1(1,:));
    post_prop_test1(2,file) = nanmean(pr_pl1(2,:));
    
    post_test_transitions(1,file) = mean(trans_short);
    post_test_transitions(2,file) = mean(trans_long);
end

figure
hold on
%x100 for % /10 since 10 trials per set of long/short,1000 is for sampling rate
dofill(1:5:5000,100/10*pre_prop_test_time_short/1000,'blue',1,60) 
dofill(1:5:5000,100/10*pre_prop_test_time_long/1000,'red',1,60)
dofill(1:5:5000,100/10*post_prop_test_time_short/1000,'green',1,60) 
dofill(1:5:5000,100/10*post_prop_test_time_long/1000,'magenta',1,60)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
title(['Novelty Preference Over Time: ' pre_files{1}(1:2)])
legend({'Pre-Short','Pre-Long','Post-Short','Post-Long'})
%%
figure
subplot(2,2,1)
bar(100*[mean(pre_prop_test,2) mean(post_prop_test,2)])
hold on
errorb(100*[mean(pre_prop_test,2) mean(post_prop_test,2)],...
    100*[std(pre_prop_test,[],2)./sqrt(length(pre_files))...
    std(post_prop_test,[],2)./sqrt(length(post_files))])
[~,p] = ttest2(pre_prop_test(1,:),post_prop_test(1,:));
if p < 0.05
    plot(1,mean(pre_prop_test(1,:))+5,'k*')
end
[~,p] = ttest2(pre_prop_test(2,:),post_prop_test(2,:));
if p < 0.05
    plot(2,mean(pre_prop_test(1,:))+5,'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
ylabel('% of time spent looking at Novel Stimulus')
yl = ylim;
ylim([50 yl(2)])
legend('Pre','Post')
title('Novelty Preference')

subplot(2,2,2)
bar([mean(pre_test_transitions,2) mean(post_test_transitions,2)])
hold on
errorb([mean(pre_test_transitions,2) mean(post_test_transitions,2)],...
    [std(pre_test_transitions,[],2)./sqrt(length(pre_files))...
    std(post_test_transitions,[],2)./sqrt(length(post_files))])
[~,p] = ttest2(pre_test_transitions(1,:),post_test_transitions(1,:));
if p < 0.05
   plot(1,mean(pre_test_transitions(1,:))+0.5,'k*');
end
[~,p] = ttest2(pre_test_transitions(2,:),post_test_transitions(2,:));
if p < 0.05
   plot(2,mean(pre_test_transitions(2,:))+0.25,'k*');
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
ylabel('# of Transitions in and out of Novel Stimulus')
legend('Pre','Post')
title('# of Transitions')

subplot(2,2,3)
bar(100*[mean(pre_prop_test1,2) mean(post_prop_test1,2)])
hold on
errorb(100*[mean(pre_prop_test1,2) mean(post_prop_test1,2)],...
    100*[std(pre_prop_test1,[],2)./sqrt(length(pre_files))...
    std(post_prop_test1,[],2)./sqrt(length(post_files))])
[~,p] = ttest2(pre_prop_test1(1,:),post_prop_test1(1,:));
if p < 0.05
    plot(1,mean(pre_prop_test1(1,:))+5,'k*')
end
[~,p] = ttest2(pre_prop_test1(2,:),post_prop_test1(2,:));
if p < 0.05
    plot(2,mean(pre_prop_test1(1,:))+5,'k*')
end
hold off
set(gca,'Xtick',[1 2])
set(gca,'XtickLabel',{'Short','Long'})
xlabel('Delay')
ylabel('% of time spent looking at Novel Stimulus')
yl = ylim;
ylim([50 yl(2)])
legend('Pre','Post')
title('Novelty Preference from 250-1000 ms')

subtitle(['Average Novelty Preference: ' pre_files{1}(1:2)])

