% analysis for alternative VPC task with no delay and long delays
% written by Seth Konig 7/2/15. Also for the RM project
%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';

% datafiles = {'PW150710.2','PW150713.2','PW150714.2'};
% sets = [7 9 11];

% datafiles = {'RR150710.2','RR150713.2','RR150714.2'};
% sets = [7 9 11];

% datafiles = {'TO150710.2','TO150713.2','TO150714.2'};
% sets = [8 10 12];

% datafiles = {'TT150713.2','TT150714.2'};
% sets = [10 12];

datafiles = {'PW150710.2','PW150713.2','PW150714.2',...
             'RR150710.2','RR150713.2','RR150714.2',...
             'TO150710.2','TO150713.2','TO150714.2',...
             'TT150713.2','TT150714.2'};
sets = [7 9 11 7 9 11 8 10 12 10 12];

% for file = 1:length(datafiles)
%     getVPC_SCW_EyeData(datafiles{file},sets(file))
% end

propleft_propright = cell(1,length(datafiles));
novel_transitions = NaN(2,length(datafiles));
for file = 1:length(datafiles)
    pl_pr = NaN(2,40);
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    trans1 = [];
    trans2 = [];
    for trial = 1:3:length(eyedat);
        x = eyedat{trial}(1,51:end);
        y = eyedat{trial}(2,51:end);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        pl_pr(1,trial) = sum(x < -7.5);
        pl_pr(2,trial) = sum(x > 7.5);
        
        timesin = findgaps(find(x < -7.5));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans1 = [trans1 count];
        
        timesin = findgaps(find(x > 7.5));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans2 = [trans2 count];
    end
    propleft_propright{file} = pl_pr; 
    novel_transitions(1,file) = mean(trans1);
    novel_transitions(2,file) = mean(trans2);
end
%%
% for familiarization phase
mean_pr_pl = NaN(2,length(datafiles));
std_pr_pl = NaN(2,length(datafiles));
for file = 1:length(datafiles)
    mean_pr_pl(1,file) = nanmean(propleft_propright{file}(1,1:3:end)./sum(propleft_propright{file}(:,1:3:end),1));
    std_pr_pl(1,file) = nanstd(propleft_propright{file}(1,1:3:end)./sum(propleft_propright{file}(:,1:3:end),1));
    mean_pr_pl(2,file) = nanmean(propleft_propright{file}(2,1:3:end)./sum(propleft_propright{file}(:,1:3:end),1));
    std_pr_pl(2,file) = nanstd(propleft_propright{file}(2,1:3:end)./sum(propleft_propright{file}(:,1:3:end),1));
end

figure
hold on
errorbar(100*mean_pr_pl(1,:),100*std_pr_pl(1,:)/sqrt(20))
errorbar(100*mean_pr_pl(2,:),100*std_pr_pl(2,:)/sqrt(20),'r')
plot([0 5],[50 50],'--k')
hold off
xlim([0 5.5])

correction = 2*mean_pr_pl./[sum(mean_pr_pl); sum(mean_pr_pl)];

% check left vs right for novel image
fam_pr_pl = NaN(2,length(datafiles));
test_rep_pr_pl1 = NaN(2,length(datafiles));
test_rep_pr_pl2 = NaN(2,length(datafiles));
for file = 1:length(datafiles)
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    fam_pr_pl(1,file) = length(find(img_pos(1,1:3:end) == -5));
    fam_pr_pl(2,file) = length(find(img_pos(1,1:3:end) == 5));
    test_rep_pr_pl1(1,file) =length(find(img_pos(1,2:3:end) == -5));
    test_rep_pr_pl1(2,file) =length(find(img_pos(1,2:3:end) == 5));
    test_rep_pr_pl2(1,file) =length(find(img_pos(1,3:3:end) == -5));
    test_rep_pr_pl2(2,file) =length(find(img_pos(1,3:3:end) == 5));
end
%%
prob_l = test_rep_pr_pl1(1,:);
prob_r = test_rep_pr_pl1(2,:);
figure
hold on
bar([nanmean(prob_l') nanmean(prob_r')]/20*100);
errorb([nanmean(prob_l') nanmean(prob_r')]/20*100,...
    [nanstd(prob_l') nanstd(prob_r')]/20*100/sqrt(5))
plot([0 3],[50 50],'k--')
hold off
set(gca,'Xtick',1:2);
set(gca,'XtickLabel',{'Left','Right'});
ylabel('% of Time on Side for Tset Phase')
ylim([40 60])
title(['Presentation Bias of Novel Test Images: ' datafiles{1}(1:2)])
%%
prop_novel1 = [];
prop_novel_corrected1= [];
prop_novel2 = [];
prop_novel_corrected2= [];
prop_novel_time1_short = zeros(length(datafiles),1000);
prop_novel_time1_long = zeros(length(datafiles),1000);
prop_novel_time2_short = zeros(length(datafiles),1000);
prop_novel_time2_long = zeros(length(datafiles),1000);

test_transitions = NaN(1,length(datafiles));

total1 = 0;
total2 = 0;
for file = 1:length(datafiles)
    pn = zeros(1,20);
    pnr = zeros(1,20);
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    if strcmpi('PW150701.2',datafiles{file})
        eyedat = [eyedat(1:37) {[NaN;NaN]}  eyedat(38:end)];
    elseif strcmpi('TT150701.2',datafiles{file})
        eyedat = [eyedat(1:38) {[NaN;NaN]}  eyedat(39:end)];
    elseif strcmpi('TT150702.2',datafiles{file})
        eyedat = [eyedat(1:38) {[NaN;NaN]}  eyedat(39:end)];
    end
    
    trans = [];
    
    for trial = 2:3:length(eyedat);
        if isnan(eyedat{trial})
            pn((trial+1)/3) = NaN;
            pnr((trial+1)/3) = NaN;
            continue
        end
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
        
        
        if img_pos(2,trial) == -12
            pn((trial+1)/3) = sum(x < -7.5)./(sum(x < -7.5) + sum(x > 7.5));
            pnr((trial+1)/3) = sum(x  < -7.5)*correction(1,file)./(sum(x < -7.5) + sum(x > 7.5));
            if trial <= 45 && trial >= 18
                prop_novel_time1_short(file,x < -7.5) =  prop_novel_time1_short(file,x < -7.5)+1;
            else
                prop_novel_time1_long(file,x < -7.5) =  prop_novel_time1_long(file,x < -7.5)+1;
            end
            total1 = total1+1;
            
            timesin = findgaps(find(x < -7.5));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            trans = [trans count];
            
        else
            pn((trial+1)/3) = sum(x > 7.5)./(sum(x < -7.5) + sum(x > 7.5));
            pnr((trial+1)/3) = sum(x > 7.5)*correction(2,file)./(sum(x < 7.5) + sum(x > 7.5));
            if trial <= 45 && trial >= 18
                prop_novel_time1_short(file,x > 7.5) =  prop_novel_time1_short(file,x > 7.5)+1;
            else
                prop_novel_time1_long(file,x > 7.5) =  prop_novel_time1_long(file,x > 7.5)+1;
            end
            total1 = total1 +1;
            
                       
            timesin = findgaps(find(x > 7.5));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            trans = [trans count];
            
        end
    end
    prop_novel1 = [prop_novel1; pn];
    prop_novel_corrected1 = [prop_novel_corrected1; pnr];
    
    for trial = 3:3:length(eyedat);
        if isnan(eyedat{trial})
            pn(trial/3) = NaN;
            pnr(trial/3) = NaN;
            continue
        end
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        
        badx = find(x < -16.5 | x > 16.5);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
        
        if img_pos(2,trial) == -12
            pn(trial/3) = sum(x < -7.5)./(sum(x < -7.5) + sum(x > 7.5));
            pnr(trial/3) = sum(x  < -7.5)*correction(1,file)./(sum(x < -7.5) + sum(x > 7.5));
            if trial <= 45 && trial >= 18
                prop_novel_time2_short(file,x < -7.5) =  prop_novel_time2_short(file,x < -7.5)+1;
            else
                prop_novel_time2_long(file,x < -7.5) =  prop_novel_time2_long(file,x < -7.5)+1;
            end
            total2 = total2+1;
            
                       
            timesin = findgaps(find(x < -7.5));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            trans = [trans count];
            
        else
            pn(trial/3) = sum(x > 7.5)./(sum(x < -7.5) + sum(x > 7.5));
            pnr(trial/3) = sum(x > 7.5)*correction(2,file)./(sum(x < 7.5) + sum(x > 7.5));
            if trial <= 45 && trial >= 18
                prop_novel_time2_short(file,x > 7.5) =  prop_novel_time2_short(file,x > 7.5)+1;
            else
                prop_novel_time2_long(file,x > 7.5) =  prop_novel_time2_long(file,x > 7.5)+1;
            end
            total2 = total2 +1;
            
                       
            timesin = findgaps(find(x > 7.5));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            trans = [trans count];
            
        end
    end
    prop_novel2 = [prop_novel2; pn];
    prop_novel_corrected2 = [prop_novel_corrected2; pnr];
    test_transitions(file) = mean(trans);
end

long_nov1 = [prop_novel1(:,1:5) prop_novel1(:,16:20)];
short_nov1 = prop_novel1(:,6:15);
long_nov_corrected1 = [prop_novel_corrected1(:,1:5) prop_novel_corrected1(:,16:20)];
short_nov_corrected1 = prop_novel_corrected1(:,6:15);

long_nov2 = [prop_novel2(:,1:5) prop_novel2(:,16:20)];
short_nov2 = prop_novel2(:,6:15);
long_nov_corrected2 = [prop_novel_corrected2(:,1:5) prop_novel_corrected2(:,16:20)];
short_nov_corrected2 = prop_novel_corrected2(:,6:15);

long_nov_by_sess1 = nanmean(long_nov1,2);
short_nov_by_sess1 =nanmean(short_nov1,2);

long_nov_by_sess_corrected1 = nanmean(long_nov_corrected1,2);
short_nov_by_sess_corrected1 = nanmean(short_nov_corrected1,2);

long_nov_by_sess2 = nanmean(long_nov2,2);
short_nov_by_sess2 = nanmean(short_nov2,2);

long_nov_by_sess_corrected2 = nanmean(long_nov_corrected2,2);
short_nov_by_sess_corrected2 = nanmean(short_nov_corrected2,2);
%%
figure
hold on
errorb([1 2],100*[nanmean(short_nov_by_sess1) nanmean(long_nov_by_sess1)],...
    100*[nanstd(short_nov_by_sess1) nanstd(long_nov_by_sess1)]./sqrt(sum(~isnan(short_nov_by_sess1))));
e1 = gca;
errorb([1 2],100*[nanmean(short_nov_by_sess_corrected1) nanmean(long_nov_by_sess_corrected1)],...
    100*[nanstd(short_nov_by_sess_corrected1) nanstd(long_nov_by_sess_corrected1)]...
    ./sqrt(sum(~isnan(short_nov_by_sess_corrected1))),'color','r');
e2 = gca;
errorb([1.5 2.5],100*[nanmean(short_nov_by_sess2) nanmean(long_nov_by_sess2)],...
    100*[nanstd(short_nov_by_sess2) nanstd(long_nov_by_sess2)]./sqrt(sum(~isnan(short_nov_by_sess2))),...
    'color','b');
e3 = gca;
errorb([1.5 2.5],100*[nanmean(short_nov_by_sess_corrected2) nanmean(long_nov_by_sess_corrected2)],...
    100*[nanstd(short_nov_by_sess_corrected2) nanstd(long_nov_by_sess_corrected2)]...
    ./sqrt(sum(~isnan(short_nov_by_sess_corrected1))),'color','m');
e4 = gca;
cp = plot([0 3],[50 50],'k--');
hold off
legend([cp e1 e2 e3 e4],{'Chance','Observed1','Bias Corrected1','Observed2','Bias Corrected2'})
set(e4,'Xtick',[1:2]);
set(e4,'XtickLabel',{'Short (~1 sec)','Long (~2 mins)'});
xlabel('Delay')
ylim([40 75])
ylabel('% of Time on Novel Stimulus')
title(['Novelty Preference: ' datafiles{1}(1:2) ' n = ' num2str(length(short_nov_by_sess2))])

%%

set_short_nov1 = NaN(1,6);
set_short_nov2 = NaN(1,6);
set_long_nov1 = NaN(1,6);
set_long_nov2 = NaN(1,6);
for s = 1:6
    ind = find(sets == s);
    set_short_nov1(s) = mean(short_nov_by_sess_corrected1(ind));
    set_short_nov2(s) = mean(short_nov_by_sess_corrected2(ind));
    set_long_nov1(s) = mean(long_nov_by_sess_corrected1(ind));
    set_long_nov2(s) = mean(long_nov_by_sess_corrected2(ind));
end

figure
subplot(2,2,1)
bar(100*set_short_nov1);
xlabel('Set#')
ylabel('% Time on Novel Stimulus')
ylim([40 75])
title('Short 1')

subplot(2,2,2)
bar(100*set_short_nov2);
xlabel('Set#')
ylabel('% Time on Novel Stimulus')
ylim([40 75])
title('Short 2')

subplot(2,2,3)
bar(100*set_long_nov1);
xlabel('Set#')
ylabel('% Time on Novel Stimulus')
ylim([40 75])
title('Long 1')

subplot(2,2,4)
bar(100*set_long_nov2);
xlabel('Set#')
ylabel('% Time on Novel Stimulus')
ylim([40 75])
title('Long 2')
subtitle('Preference across all monkeys and sets (n = 1-3)')
%%
figure
hold on
dofill(1:5:5000,prop_novel_time1_short/(5*(total1/length(datafiles))),'blue',1,50)
dofill(1:5:5000,prop_novel_time1_long/(5*(total1/length(datafiles))),'green',1,50)
dofill(1:5:5000,prop_novel_time2_short/(5*(total2/length(datafiles))),'red',1,50)
dofill(1:5:5000,prop_novel_time2_long/(5*(total2/length(datafiles))),'magenta',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
legend('Test 1 Short','Test 1 Long','Test 2 Short','Test 2 Long')
title(['Novelty Preference Over Time: ' datafiles{1}(1:2)])
%%
figure
hold on
bar([mean(novel_transitions') mean(test_transitions)])
errorb([mean(novel_transitions') mean(test_transitions)],...
    [std(novel_transitions') std(test_transitions)]./sqrt(length(test_transitions)))
hold off
set(gca,'Xtick',1:3)
set(gca,'XtickLabel',{'Familarization Left','Familarization Right','Test'})
ylabel('# of Transitions')
title('VPC SC Wide')
%%
pw_short = prop_novel_time1_short(1:3,:)+prop_novel_time2_short(1:3,:);
pw_long = prop_novel_time1_long(1:3,:)+prop_novel_time2_long(1:3,:);
rr_short = prop_novel_time1_short(4:6,:)+prop_novel_time2_short(4:6,:);
rr_long = prop_novel_time1_long(4:6,:)+prop_novel_time2_long(4:6,:);
to_short = prop_novel_time1_short(7:9,:)+prop_novel_time2_short(7:9,:);
to_long = prop_novel_time1_long(7:9,:)+prop_novel_time2_long(7:9,:);
tt_short = prop_novel_time1_short(10:11,:)+prop_novel_time2_short(10:11,:);
tt_long = prop_novel_time1_long(10:11,:)+prop_novel_time2_long(10:11,:);
%%
figure
subplot(1,2,1)
hold on
dofill(1:5:5000,pw_short/(10*(total1/length(datafiles))),'blue',1,50)
dofill(1:5:5000,rr_short/(10*(total1/length(datafiles))),'red',1,50)
dofill(1:5:5000,to_short/(10*(total1/length(datafiles))),'green',1,50)
dofill(1:5:5000,tt_short/(10*(total1/length(datafiles))),'magenta',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
legend('Vivian','Red','Tobii','Timmy')
title('Short')

subplot(1,2,2)
hold on
dofill(1:5:5000,pw_long/(10*(total1/length(datafiles))),'blue',1,50)
dofill(1:5:5000,rr_long/(10*(total1/length(datafiles))),'red',1,50)
dofill(1:5:5000,to_long/(10*(total1/length(datafiles))),'green',1,50)
dofill(1:5:5000,tt_long/(10*(total1/length(datafiles))),'magenta',1,50)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
legend('Vivian','Red','Tobii','Timmy')
title('Long')