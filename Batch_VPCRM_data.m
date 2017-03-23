% Run analysis for all VPC data for Relational Memory (RM) project
% written by Seth Konig June, 2015. Note by Seth Konig 6/9/16 probably
% bugy definitley not written efficiently. See VPCRM_Pre_vs_Post.

%set 24 images 53 and 54 are the same

%---Vivian Pre-Lesion---%
% datafiles = {'PW150609.2','PW150610.2','PW150611.2','PW150612.2',...
%              'PW150615.2','PW150616.2','PW150617.2','PW150618.2',...
%              'PW150619.2','PW150622.2','PW150623.2','PW150624.2',...
%              'PW150625.2','PW150626.2','PW150629.2','PW150630.2'};
% setnums = [41 42 43 44 ...
%            45 46 47 48 ...
%            49 50 51 52 ...
%            53 54 55 56];
       
%---Vivian Post-lesion---%
% datafiles = {'PW160425.2','PW160426.2','PW160427.2','PW160428.2','PW160429.2',...
%              'PW160502.2','PW160503.2','PW160504.2','PW160505.2','PW160506.2',...
%              'PW160509.2','PW160513.2','PW160516.2','PW160518.2','PW160519.2'};
% setnums = [20 21 22 23 24,...
%            25 26 27 28 29,...
%            30 32 33 34 35];


%---Manfred Pre-lesion---%
datafiles = {'MF170320.3','MF170321.2'};
setnums = [41 42];

%---Red Pre-Lesion---%
% datafiles = {'RR150609.2','RR150610.2','RR150611.2','RR150612.2',...
%              'RR150615.2','RR150616.2','RR150617.2','RR150618.2',...
%              'RR150619.2','RR150622.2','RR150623.2','RR150624.2',...
%              'RR150625.2','RR150626.2','RR150629.2','RR150630.2'};
% setnums = [20 21 22 23 ...
%            24 25 26 27 ...
%            28 29 30 31 ...
%            32 33 34 35];


%---Red Post-Lesion---%

% datafiles = {'RR160523.2','RR160526.2','RR160531.2','RR160601.2',...
%              'RR160602.2','RR160603.2','RR160606.2','RR160607.2',...
%              'RR160608.2','RR160609.2','RR160610.2','RR160613.2',...
%              'RR160614.2','RR160615.2'};
% 
% setnums = [41 42 43 44 ...
%            45 46 47 48 ...
%            49 50 51 52 ...
%            57 58];

% datafiles = {'TT150609.2','TT150610.2','TT150611.2','TT150612.2',...
%              'TT150615.2','TT150616.2','TT150617.2','TT150618.2'....
%              'TT150622.2','TT150623.2','TT150624.2','TT150625.2',...
%              'TT150626.2','TT150629.2','TT150630.2'};
% setnums = [20 21 22 23 ...
%            24 25 26 27 ...
%            28 29 30 31 ...
%            32 33 34];
%
% datafiles = {'TO150609.2','TO150610.2','TO150611.2','TO150612.2',...
%     'TO150615.2','TO150616.2','TO150617.2','TO150618.2'....
%     'TO150622.2','TO150623.2','TO150624.2','TO150625.2',...
%     'TO150626.2','TO150629.2','TO150630.2'};
% setnums = [20 21 22 23 ...
%     24 25 26 27 ...
%     28 29 30 31 ...
%     32 33 34];

% datafiles = {'TO150609.2','TO150610.2','TO150611.2','TO150612.2',...
%     'TO150615.2','TO150616.2','TO150617.2','TO150618.2'....
%     'TO150622.2','TO150623.2','TO150624.2','TO150625.2',...
%     'TO150626.2','TO150629.2','TO150630.2',...
%     'RR150609.2','RR150610.2','RR150611.2','RR150612.2',...
%     'RR150615.2','RR150616.2','RR150617.2','RR150618.2',...
%     'RR150619.2','RR150622.2','RR150623.2','RR150624.2',...
%     'RR150625.2','RR150626.2','RR150629.2','RR150630.2',...
%     'PW150609.2','PW150610.2','PW150611.2','PW150612.2',...
%     'PW150615.2','PW150616.2','PW150617.2','PW150618.2',...
%     'PW150619.2','PW150622.2','PW150623.2','PW150624.2',...
%     'PW150625.2','PW150626.2','PW150629.2','PW150630.2'};
% setnums = [20 21 22 23 ...
%     24 25 26 27 ...
%     28 29 30 31 ...
%     32 33 34 ...
%     20 21 22 23 ...
%     24 25 26 27 ...
%     28 29 30 31 ...
%     32 33 34 35 ...
%     41 42 43 44 ...
%     45 46 47 48 ...
%     49 50 51 52 ...
%     53 54 55 56];

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\';

for file = 1:length(datafiles)
    getVPC_EyeData(datafiles{file},setnums(file));
end
%%
novel_transitions = NaN(2,length(datafiles));
propleft_propright = cell(1,length(datafiles));
all_pupil = [];
for file = 1:length(datafiles)
    pl_pr = NaN(2,40);
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    trans1 = [];
    trans2 = [];
    for trial = 1:length(eyedat);
        x = eyedat{trial}(1,:);
        y = eyedat{trial}(2,:);
        img_on = round(per(trial).alltim(per(trial).allval == 23)/5);
        all_pupil = [all_pupil pupildata{trial}(img_on-100:img_on+500)];
        
        
        badx = find(x < -8.75 | x > 8.75);
        bady = find(y < -4.5 | y > 4.5);
        
        badind = union(badx,bady);
        x(badind) = [];
        y(badind) = [];
        
        pl_pr(1,trial) = sum(x < 0);
        pl_pr(2,trial) = sum(x > 0);
        
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
end
%%
% for familiarization phase
mean_pr_pl = NaN(2,length(datafiles));
std_pr_pl = NaN(2,length(datafiles));
for file = 1:length(datafiles)
    mean_pr_pl(1,file) = mean(propleft_propright{file}(1,1:2:end));
    std_pr_pl(1,file) = std(propleft_propright{file}(1,1:2:end));
    mean_pr_pl(2,file) = mean(propleft_propright{file}(2,1:2:end));
    std_pr_pl(2,file) = std(propleft_propright{file}(2,1:2:end));
end

figure
hold on
errorbar(100*mean_pr_pl(1,:)/1000,100*std_pr_pl(1,:)/1000/sqrt(20))
errorbar(100*mean_pr_pl(2,:)/1000,100*std_pr_pl(2,:)/1000/sqrt(20),'r')
plot([0 5],[50 50],'--k')
hold off
xlim([0.5 9.5])

figure
hold on
bar(nanmean(mean_pr_pl')/10)
errorb(nanmean(mean_pr_pl')/10,nanstd(mean_pr_pl')/sqrt(5)/10)
plot([0 3],[50 50],'k--')
hold off
set(gca,'XTick',[1:2])
set(gca,'XtickLabel',{'Left','Right'})
ylabel('% of Time on Side for Familarization Phase')
ylim([35 65])
title(['Left/Right Bias: ' datafiles{1}(1:2)])
%%
correction = 2*mean_pr_pl./[sum(mean_pr_pl); sum(mean_pr_pl)];

% check left vs right for novel image
fam_pr_pl = NaN(2,length(datafiles));
test_rep_pr_pl = NaN(2,length(datafiles));
test_nov_pr_pl = NaN(2,length(datafiles));
for file = 1:length(datafiles)
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    fam_pr_pl(1,file) = length(find(img_pos(1,1:2:end) == -5));
    fam_pr_pl(2,file) = length(find(img_pos(1,1:2:end) == 5));
    test_rep_pr_pl(1,file) =length(find(img_pos(1,2:2:end) == -5));
    test_rep_pr_pl(1,file) =length(find(img_pos(1,2:2:end) == 5));
    test_nov_pr_pl(1,file) =length(find(img_pos(2,2:2:end) == -5));
    test_nov_pr_pl(1,file) =length(find(img_pos(2,2:2:end) == 5));
end

prob_l = test_rep_pr_pl(1,:);
prob_r = 20-test_rep_pr_pl(1,:);
figure
hold on
bar([nanmean(prob_l') nanmean(prob_r')]/20*100);
errorb([nanmean(prob_l') nanmean(prob_r')]/20*100,...
    [nanstd(prob_l') nanstd(prob_r')]/20*100/sqrt(5))
plot([0 3],[50 50],'k--')
hold off
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Left','Right'})
ylabel('% of Time on Side for Tset Phase')
ylim([40 60])
title(['Presentation Bias of Novel Test Images: ' datafiles{1}(1:2)])
%%
prop_novel = [];
prop_novel_corrected = [];
prop_novel_time_short = zeros(length(datafiles),1000);
prop_novel_time_long = zeros(length(datafiles),1000);

test_transitions = NaN(1,length(datafiles));

for file = 1:length(datafiles)
    pn = zeros(1,20);
    pnr = zeros(1,20);
    load([data_dir datafiles{file}(1:8) '_' datafiles{file}(10) '-fixation.mat']);
    
    trans = [];
    for trial = 2:2:length(eyedat);
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
        elseif length(x) < 1000 %sometimes cortex isn't quite perfect at timing
            x = [x zeros(1,1000-length(x))];
            y = [y zeros(1,1000-length(y))];
        end
        
        if img_pos(2,trial) == -5
            pn(trial/2) = sum(x < 0);
            pnr(trial/2) = sum(x < 0)*correction(1,file);
            if trial <= 30 && trial >= 12
                prop_novel_time_short(file,x <0) =  prop_novel_time_short(file,x <0)+1;
                
            else
                prop_novel_time_long(file,x <0) =  prop_novel_time_long(file,x <0)+1;
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
            pn(trial/2) = sum(x > 0);
            pnr(trial/2) = sum(x > 0)*correction(2,file);
            
            if trial <= 30 && trial >= 12
                prop_novel_time_short(file,x > 0 ) =  prop_novel_time_short(file,x >0)+1;
            else
                prop_novel_time_long(file,x > 0) =  prop_novel_time_long(file,x >0)+1;
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
    prop_novel = [prop_novel; pn];
    prop_novel_corrected = [prop_novel_corrected; pnr];
    test_transitions(file) = mean(trans);
end
long_nov = [prop_novel(:,1:5) prop_novel(:,16:20)];
short_nov = prop_novel(:,6:15);
long_nov_corrected = [prop_novel_corrected(:,1:5) prop_novel_corrected(:,16:20)];
short_nov_corrected = prop_novel_corrected(:,6:15);

long_nov_by_sess = mean(long_nov,2)/10;
short_nov_by_sess = mean(short_nov,2)/10;

long_nov_by_sess_corrected = mean(long_nov_corrected,2)/10;
short_nov_by_sess_corrected = mean(short_nov_corrected,2)/10;

% n_long = sampsizepwr('z',[mean(long_nov_by_sess),std(long_nov_by_sess)],500);
% n_short= sampsizepwr('z',[mean(short_nov_by_sess),std(short_nov_by_sess)],500);

figure
hold on
cp = plot([0 3],[50 50],'k--');
plot([1 2],[nanmean(short_nov_by_sess) nanmean(long_nov_by_sess)],'k.')
plot([1 2],[nanmean(short_nov_by_sess_corrected) nanmean(long_nov_by_sess_corrected)],'r.')
errorb([nanmean(short_nov_by_sess) nanmean(long_nov_by_sess)],...
    [nanstd(short_nov_by_sess) nanstd(long_nov_by_sess)]./sqrt(sum(~isnan(short_nov_by_sess))));
errorb([nanmean(short_nov_by_sess_corrected) nanmean(long_nov_by_sess_corrected)],...
    [nanstd(short_nov_by_sess_corrected) nanstd(long_nov_by_sess_corrected)]...
    ./sqrt(sum(~isnan(short_nov_by_sess_corrected))),'color','r');
hold off
legend({'Chance','Observed','Bias Corrected'})
% legend('Chance','Observed','Bias Corrected')
set(gca,'Xtick',1:2)
set(gca,'XtickLabel',{'Short (~10 secs)','Long (~2 mins)'})
xlabel('Delay')
ylim([45 80])
ylabel('% of Time on Novel Stimulus')
title(['Novelty Preference: ' datafiles{1}(1:2)])

%%
figure
hold on
dofill(1:5:5000,prop_novel_time_short/100,'blue',1,50)
dofill(1:5:5000,prop_novel_time_long/100,'red',1,50)

plot([0 5000],[50 50],'k--')
hold off
xlabel('Time (ms)')
ylabel('% time looking at novel stimulus')
title(['Novelty Preference Over Time: ' datafiles{1}(1:2)])
legend({'Short','Long'})


%%
% %% Average Left Right bias by Monkey
% pw_plr = [];
% rr_plr = [];
% to_plr = [];
% tt_plr = [];
% for copy = 1:4
%     pw_plr(:,copy) = mean(propleft_propright{4*(copy-1) + 1},2)/10;
%     rr_plr(:,copy)= mean(propleft_propright{4*(copy-1) + 2},2)/10;
%     to_plr(:,copy) = mean(propleft_propright{4*(copy-1) + 3},2)/10;
%     tt_plr(:,copy) = mean(propleft_propright{4*(copy-1) + 4},2)/10;
% end
%
% figure
% hold on
% errorbar(1,mean(pw_plr(1,:)),std(pw_plr(1,:))'./sqrt(size(pw_plr,2)),'b')
% errorbar(1,mean(pw_plr(2,:)),std(pw_plr(2,:))'./sqrt(size(pw_plr,2)),'r')
% errorbar(2,mean(rr_plr(1,:)),std(rr_plr(1,:))'./sqrt(size(rr_plr,2)),'b')
% errorbar(2,mean(rr_plr(2,:)),std(rr_plr(2,:))'./sqrt(size(rr_plr,2)),'r')
% errorbar(3,mean(to_plr(1,:)),std(to_plr(1,:))'./sqrt(size(to_plr,2)),'b')
% errorbar(3,mean(to_plr(2,:)),std(to_plr(2,:))'./sqrt(size(to_plr,2)),'r')
% errorbar(4,mean(tt_plr(1,:)),std(tt_plr(1,:))'./sqrt(size(tt_plr,2)),'b')
% errorbar(4,mean(tt_plr(2,:)),std(tt_plr(2,:))'./sqrt(size(tt_plr,2)),'r')
% hold off
% set(gca,'Xtick',1:4)
% set(gca,'XtickLabel',{'PW','RR','TO','TT'})
% xlabel('Monkey')
% ylabel('% of Time')
% legend('Left','Right')
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
title('VPC Flickr')