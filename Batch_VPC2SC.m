% analysis for alternative VPC task with no delay and long delays
% written by Seth Konig 7/2/15. Also for the RM project
%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';
%
%---Pre-Lesion---%
% datafiles = {'PW150701.2','PW150706.2','PW150707.2'};
% sets =  [1 2 5];

% datafiles = {'RR150701.2','RR150706.2','RR150707.2'};
% sets = [1 2 5];

% datafiles = {'TO150702.2','TO150706.2','TO150707.2'};
% sets =  [3 4 6];

% datafiles = {'MF170412.2','MF170413.3'};
% sets = [3 4];

%for across sets and monkeys
% datafiles = {'PW150701.2','PW150706.2','PW150707.2',...
%     'RR150701.2','RR150706.2','RR150707.2',...
%     'TO150702.2','TO150706.2','TO150707.2',...
%     'MF170412.2','MF170413.3'};
% sets =  [1 2 5 ...
%     1 2 5 ...
%     3 4 6 ...
%     3 4];

%---Post-Lesion---%
datafiles = {'PW160520.1','PW160524.2','PW160525.2'}; %task may not have worked right
sets =  [3 4 6];


%  datafiles = {'RR160616.2','RR160617.2','RR160620.2'};
%  sets = [3 4 6];

% datafiles = {'TO170810.2','TO170811.2','TO170814.2'};
% sets =  [1 2 5];

%---Timmy---%
% datafiles = {'TT150701.2','TT150702.2','TT150706.2','TT150707.2'};
% sets = [1 3 2 6];




for file = 1:length(datafiles)
    getVPC_SC_EyeData(datafiles{file},sets(file))
end
%%
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
        
        badx = find(x < -8.75 | x > 8.75);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        pl_pr(1,trial) = sum(x < 0);
        pl_pr(2,trial) = sum(x > -0);
        
        timesin = findgaps(find(x < -0.5));
        count = 0;
        for t = 1:size(timesin,1)
            ind = timesin(t,:) ~= 0;
            if sum(ind) > 10
                count = count+1;
            end
        end
        trans1 = [trans1 count];
        
        timesin = findgaps(find(x > 0.5));
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
    
    %error in block sturcture resulted in skipping a trial when break
    %fixation occured in block 45
    
    if strcmpi('PW150701.2',datafiles{file})
        propleft_propright{file} = [ propleft_propright{file}(:,1:37) [NaN;NaN]  propleft_propright{file}(:,38:end)];
    elseif strcmpi('TT150701.2',datafiles{file})
        propleft_propright{file} = [ propleft_propright{file}(:,1:38) [NaN;NaN]  propleft_propright{file}(:,39:end)];
    elseif strcmpi('TT150702.2',datafiles{file})
        propleft_propright{file} = [ propleft_propright{file}(:,1:38) [NaN;NaN]  propleft_propright{file}(:,39:end)];
    end
end

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
% set(gca,'Xtick',1:2);
% set(gca,'XtickLabel',{'Left','Right'});
ylabel('% of Time on Side for Tset Phase')
ylim([40 60])
title(['Presentation Bias of Novel Test Images: ' datafiles{1}(1:2)])
%%
prop_novel1 = [];
prop_novel_corrected1= [];
prop_novel2 = [];
prop_novel_corrected2= [];
prop_novel_time_short1 = zeros(length(datafiles),1000);
prop_novel_time_short2 = zeros(length(datafiles),1000);
prop_novel_time_long1 = zeros(length(datafiles),1000);
prop_novel_time_long2 = zeros(length(datafiles),1000);
totalshort1 = 0;
totalshort2 = 0;
totallong1 = 0;
totallong2 = 0;


test_transitions = NaN(1,length(datafiles));

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
        
        badx = find(x < -8.75 | x > 8.75);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
        
        
        if img_pos(2,trial) == -5
            pn((trial+1)/3) = sum(x < 0)./(sum(x < -0) + sum(x > 0));
            pnr((trial+1)/3) = sum(x  < 0)*correction(1,file)./(sum(x < -0) + sum(x > 0));
            
            if trial <= 15 || trial > 45
                prop_novel_time_long1(file,x < 0) =  prop_novel_time_long1(file,x < 0)+1;
                totallong1 = totallong1+1;
            else
                prop_novel_time_short1(file,x < 0) =  prop_novel_time_short1(file,x < 0)+1;
                totalshort1 = totalshort1+1;
            end
            
            
            timesin = findgaps(find(x< -0.5));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            trans = [trans count];
        else
            pn((trial+1)/3) = sum(x > 0)./(sum(x < -0) + sum(x > 0));
            pnr((trial+1)/3) = sum(x > 0)*correction(2,file)./(sum(x < 0) + sum(x > 0));
            
            if trial <= 15 || trial > 45
                prop_novel_time_long1(file,x < 0) =  prop_novel_time_long1(file,x < 0)+1;
                totallong1 = totallong1+1;
            else
                prop_novel_time_short1(file,x < 0) =  prop_novel_time_short1(file,x < 0)+1;
                totalshort1 = totalshort1+1;
            end
            
            
            timesin = findgaps(find(x > 0.5));
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
        
        badx = find(x < -8.75 | x > 8.75);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        if length(x) > 1000
            x = x(1:1000);
            y = y(1:1000);
        end
        
        if img_pos(2,trial) == -5
            pn(trial/3) = sum(x < -0)./(sum(x < -0) + sum(x > 0));
            pnr(trial/3) = sum(x < -0)*correction(1,file)./(sum(x < -0) + sum(x > 0));
            
            if trial <= 15 || trial > 45
                prop_novel_time_long2(file,x < 0) =  prop_novel_time_long2(file,x < 0)+1;
                totallong2 = totallong2+1;
            else
                prop_novel_time_short2(file,x < 0) =  prop_novel_time_short2(file,x < 0)+1;
                totalshort2= totalshort2+1;
            end
            
            
            
            timesin = findgaps(find(x< -0.5));
            count = 0;
            for t = 1:size(timesin,1)
                ind = timesin(t,:) ~= 0;
                if sum(ind) > 10
                    count = count+1;
                end
            end
            trans = [trans count];
        else
            pn(trial/3) = sum(x > 0)./(sum(x < -0) + sum(x > 0));
            pnr(trial/3) = sum(x >0)*correction(2,file)./(sum(x < -0) + sum(x > 0));
            
            if trial <= 15 || trial > 45
                prop_novel_time_long2(file,x < 0) =  prop_novel_time_long2(file,x < 0)+1;
                totallong22 = totallong2+1;
            else
                prop_novel_time_short2(file,x < 0) =  prop_novel_time_short2(file,x < 0)+1;
                totalshort2 = totalshort2+1;
            end
            
            
            timesin = findgaps(find(x > 0.5));
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
% set(e4,'Xtick',[1:2]);
% set(e4,'XtickLabel',{'Short (~1 sec)','Long (~2 mins)'});
xlabel('Delay')
ylim([40 75])
ylabel('% of Time on Novel Stimulus')
title(['Novelty Preference: n = ' num2str(length(short_nov_by_sess2))])

%%
figure
hold on
errorb([0.75 1.75],100*[nanmean(short_nov_by_sess1) nanmean(long_nov_by_sess1)],...
    100*[nanstd(short_nov_by_sess1) nanstd(long_nov_by_sess1)]./sqrt(sum(~isnan(short_nov_by_sess1))));
errorb([1.25 2.25],100*[nanmean(short_nov_by_sess2) nanmean(long_nov_by_sess2)],...
    100*[nanstd(short_nov_by_sess2) nanstd(long_nov_by_sess2)]./sqrt(sum(~isnan(short_nov_by_sess2))),...
    'color','b');
cp = plot([0 3],[50 50],'k--');
plot(1,50,'k')
plot(1,50,'b')
hold off
legend('Chance','Observed1','Observed2')
set(gca,'Xtick',[1:2]);
set(gca,'XtickLabel',{'Short (~1 sec)','Long (~2 mins)'});
xlabel('Delay')
ylim([40 75])
ylabel('% of Time on Novel Stimulus')
title(['VPC2SC: n = ' num2str(length(short_nov_by_sess2))])
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
dofill(1:5:5000,prop_novel_time_short1./totalshort1,'black',1,25)
dofill(1:5:5000,prop_novel_time_short2./totalshort1,'green',1,25)
dofill(1:5:5000,prop_novel_time_long1./totalshort1,'blue',1,25)
dofill(1:5:5000,prop_novel_time_long2./totalshort1,'magenta',1,25)
plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
legend('Short 1','Short 2','Long 1','Long 2')
title(['VPC2SC: Novelty Preference Over Time'])
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
title('VPC 2 SC')