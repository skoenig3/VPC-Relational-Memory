function getVPC_EyeData(datafile,itmnum)
%modified from getSCM_EyeDat by Seth Konig 6/9/2015
%function imports VPC eye data and calibrates it. 

figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\VPC Figures\';
samprate = 5;

ITMFile = ['C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Item and CND Files\2im0' num2str(itmnum) '.itm'];
CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Item and CND files\2im60.cnd';

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

itmposx = zeros(1,size(itmfil,1)-8);
for i = 7:size(itmfil,1)-1
    str = textscan(itmfil(i,:),'%d');
   itmposx(i-6) =  str{1}(3); 
end

if strcmpi(datafile(1:2),'PW')
    datafile = ['Z:\Cortex Data\Vivian\' datafile];
elseif strcmpi(datafile(1:2),'TT')
    datafile = ['Z:\Cortex Data\Timmy\' datafile];
elseif strcmpi(datafile(1:2),'RR')
    datafile = ['Z:\Cortex Data\Red\' datafile];
elseif strcmpi(datafile(1:2),'TO')
    datafile = ['Z:\Cortex Data\Tobii\' datafile];
elseif strcmpi(datafile(1:2),'MF')
    datafile = ['Z:\Cortex Data\Manfred\' datafile];
end

%sub function convert raw eye tracking data into x & y coordinates
[time_arr,event_arr,eog_arr,epp_arr,~,~]  = get_ALLdata(datafile);
numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 1010)) ~=0
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
            perendind = find(event_arr(:,rptlop) == 24);
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop)-100;
            endtimdum = time_arr(perendind,rptlop);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                clrchgind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
            end
        end
    end
end
%clear cnd
numrpt = size(per,2);
for rptlop = 1:numrpt
    cnd(rptlop)=per(rptlop).cnd;
end
evnnmb=2:2:size(eog_arr,1);
oddnmb=1:2:size(eog_arr,1);
clear x y
cndlst=unique(cnd);
for k=1:length(cndlst)
    cndind=find(cnd==cndlst(k));
    allind=clrchgind(cndind);
    for l=1:length(allind)
        x{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,oddnmb),allind(l)));
        y{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,evnnmb),allind(l)));
    end
end
%remove outlying points when calculating average eye position @ each location
for k=1:length(x)
    x{k}=x{k}(find(x{k}<mean(x{k}+std(x{k})) & x{k}>mean(x{k}-std(x{k}))));
    y{k}=y{k}(find(x{k}<mean(x{k}+std(x{k})) & x{k}>mean(x{k}-std(x{k}))));
    x{k}=x{k}(find(y{k}<mean(y{k}+std(y{k})) & y{k}>mean(y{k}-std(y{k}))));
    y{k}=y{k}(find(y{k}<mean(y{k}+std(y{k})) & y{k}>mean(y{k}-std(y{k}))));
end
clear meanx meany
for k=1:length(x)
    meanx(k)=mean(x{k});
end
for k=1:length(y)
    meany(k)=mean(y{k});
end
%%

% old calibration code--removed 5/20/15
% clear x y
% x=meanx; y=meany;
% meanxorigin = mean([x(6) x(2) x(1) x(5) x(9) ],2);
% xscale = mean([6/(x(8)-meanxorigin) 3/(x(4)-meanxorigin) 3/(abs(x(3)-meanxorigin)) 6/(abs(x(7)-meanxorigin))],2);
% meanyorigin = mean([y(7) y(3) y(1) y(4) y(8) ],2);
% yscale = mean([6/(y(6)-meanyorigin) 3/(y(2)-meanyorigin) 3/(abs(y(5)-meanyorigin)) 6/(abs(y(9)-meanyorigin))],2);

%new calibration code
controlx = [0 0 -3 3 0 0 -6 6 0];
controly = [0 3 0 0 -3 6 0 0 -6];
tform = get_calibration_fcn([controlx; controly],[meanx;meany]);

%%
figure
subplot(1,2,1)
hold on
for i = 1:9
    xi = x{i};
    yi = y{i};
    plot(xi,yi,'*')
end
hold off
title('Raw Calibration Data')

subplot(1,2,2)
hold on
hold on
for i = 1:9
    xi = x{i};
    yi = y{i};
    [xi,yi] = tformfwd(tform,xi,yi);
    plot(xi,yi,'*')
end
plot(controlx,controly,'k+','markersize',5)
hold off
title('Post-Hoc Calibration')
save_and_close_fig(figure_dir,['Estimated_Calibration_'  datafile(end-9:end-2)])

%%


numrpt = size(event_arr,2);
valrptcnt = 0;
clear per vpcind

new_eog_arr=[];
new_epp_arr = [];
for rptlop = 1:numrpt
    if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
        if size(find(event_arr(:,rptlop) == 200)) ~=0
            perbegind = find(event_arr(:,rptlop) == 23,1,'first');
            perendind = find(event_arr(:,rptlop) == 24,1,'first');
            cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
            blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
            begtimdum = time_arr(perbegind,rptlop);
            endtimdum = time_arr(perendind,rptlop);
            if endtimdum > begtimdum
                valrptcnt = valrptcnt + 1;
                vpcind(valrptcnt)=rptlop;
                per(valrptcnt).begsmpind = begtimdum;
                per(valrptcnt).endsmpind = endtimdum;
                per(valrptcnt).begpos = 1;
                per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                per(valrptcnt).blk = event_arr(blknumind,rptlop);
                per(valrptcnt).allval = event_arr(:,rptlop);
                per(valrptcnt).alltim = time_arr(:,rptlop);
                new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                new_epp_arr = cat(2,new_epp_arr,epp_arr(:,rptlop));
                cnd = find(event_arr(:,rptlop) >= 1010 & event_arr(:,rptlop) < 2000);
                type = event_arr(cnd-1,rptlop);
            end
        end
    end
end

eyedat = [];
if ~isempty(epp_arr);
    pupildata = cell(1,length(per)); %horizontal pupil diameter
end
for trlop=1:size(per,2)
    trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
    horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
    vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
    picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to eye scan start
    picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to eye scan start
    
    if picend > length(horeog)*5
        picend =  length(horeog)*5;
    end
    eyedat{trlop}(1,:) = (horeog(ceil(picstart/5):floor(picend/5)));% .* xscale; %removed 5/20/15
    eyedat{trlop}(2,:) = (vrteog(ceil(picstart/5):floor(picend/5)));% .* yscale; %removed 5/20/15
    
    try
        trlepp=new_epp_arr(~isnan(new_epp_arr(:,trlop)),trlop); % epp for this trial
        pupildata{trlop} = trlepp(2:2:2*floor(picend/samprate));%odd indexes contains nothing but noise
    catch
        if picend > 20000
            pupildata{trlop} = trlepp(2:2:end); %epp for this trial
        else
            try
                pupildata{trlop} = trlepp(2:2:(2*floor(picend/samprate)-1));%odd indexes contains nothing but noise
            catch
                error('Unknown reason EOG and pic duration do not line up')
            end
        end
    end
end


%---Recalibrate and automatically scale eye data---%
% added 5/20/15
for eye = 1:length(eyedat)
    if all(~isnan(eyedat{eye}))
        x = eyedat{eye}(1,:);
        y = eyedat{eye}(2,:);
        [x,y] = tformfwd(tform,x,y);
        eyedat{eye} = [x;y];
    end
end

img_cnd = zeros(2,40);
img_trial = 1;
for trial = 1:length(per);
    cndline = textscan(cndfil(per(trial).cnd+1,:),'%d');
    imgnums = [cndline{1}(end-1) cndline{1}(end)];
    img_cnd(:,trial) = imgnums';
end

imgnum = zeros(2,40); 
for trial = 1:size(img_cnd,2)
    itmline1 = textscan(itmfil(img_cnd(1,trial)+6,:),'%s');
    imgnum1 = str2double(itmline1{1}{end}(11:end-4));
    imgnum(1,trial) = imgnum1;
    
    itmline2 = textscan(itmfil(img_cnd(2,trial)+6,:),'%s');
    imgnum2 = str2double(itmline2{1}{end}(11:end-4));
    imgnum(2,trial) = imgnum2;
end

img_pos = itmposx(img_cnd);
    

%check trials per block
% block_trials = zeros(1,79);
% for block = 1:79
%     [~,colind] = find(event_arr == 500+block);
%     subtract = 0; 
%     for c = 1:length(colind)
%         if isempty(find(event_arr(:,colind(c)) == 3)) && isempty(find(event_arr(:,colind(c)) == 200))
%             subtract = subtract + 1; 
%         end
%     end
%     block_trials(block) = length(colind)-subtract; 
% end
    

% %remove eye data outside of image
% for i = 1:size(eyedat,2);
%     x = eyedat{i}(1,:)*24+400;
%     y = eyedat{i}(2,:)*24+300;
%     %code shouldn't be necessary since image start is already time 0
%     %     tstart = find(x > 300 & x < 500 & y > 200 & y < 400);
%     %     tstart = tstart(1);
%     %     x = x(tstart:end); y = y(tstart:end);
%     if length(x) > 2000; %1st 11 seconds of data, ~10 seconds but 1 second buffer for accidental lookaways, etc
%         x(2001:end) = []; y(2001:end) = [];
%     end
%     badx = find(x < -25 | x > 825); %~1 dva leave margin of error, was 2 dva SDK 5/26 but calibration looks better now
%     x(badx) = []; y(badx) = [];
%     bady = find(y < -25 | y > 625); %~1 dva margin of error, was 2 dva SDK 5/26 but calibration looks better now
%     x(bady) = []; y(bady) = [];
%     x = (x-400)/24; y = (y-300)/24;
%     eyedat{i} = [x;y];
% end
% end

save(['C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Eye Data\'...
    datafile(end-9:end-2) '_' datafile(end) '-fixation'],'eyedat','per','imgnum',...
    'img_cnd','img_pos','pupildata');