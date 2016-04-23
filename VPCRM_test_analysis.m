datafile = 'VPCRMTES.1';


ITMFile = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Item and CND Files\2im020.itm';
CNDFile = 'C:\Users\seth.koenig\Documents\MATLAB\VPC Relational Memory\Item and CND files\2im60.cnd';

[time_arr,event_arr,eog_arr,epp_arr, header,trialcount]  = get_ALLdata(datafile);


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

itmlist = zeros(size(cndfil,1)-1,2);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1,1) = str{1}(end-1);
    itmlist(i-1,2) = str{1}(end);
end

itmposx = zeros(1,size(itmfil,1)-8);
for i = 7:size(itmfil,1)-1
    str = textscan(itmfil(i,:),'%d');
   itmposx(i-6) =  str{1}(3); 
end

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 > 1 %1st block is color change always
        if event_arr(find(event_arr(:,rptlop)>1000,1,'last'),rptlop)-1000 <= 9 %clrchng trials
            if size(find(event_arr(:,rptlop) == 3)) ~=0
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
end

numrpt = size(per,2);
for rptlop = 1:numrpt
    clrcnd(rptlop) = per(rptlop).cnd-1000;
end
count = NaN(1,9);
for item = 1:9;
   count(item) = sum(clrcnd == item);  
end
%%
numrpt = size(event_arr,2);
valrptcnt = 0;
clear per clrchgind
for rptlop = 1:numrpt
    if event_arr((find(event_arr(:,rptlop)>500,1,'last')),rptlop)-500 > 1 %1st block is color change always
        if event_arr(find(event_arr(:,rptlop)>1000,1,'last'),rptlop)-1000 > 9 %clrchng trials
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
                    
                    cross_on = find(event_arr(:,rptlop) == 9);
                    cross_off = find(event_arr(:,rptlop) == 10);
                    img_on = find(event_arr(:,rptlop) == 23);
                    img_off = find(event_arr(:,rptlop) == 24);
                    
                    per(valrptcnt).crossfixdur = time_arr(cross_off,rptlop)-time_arr(cross_on,rptlop);
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop)-1000;
                    per(valrptcnt).imgdur = time_arr(img_off,rptlop)-time_arr(img_on,rptlop);
                end
            end
        end
    end
end


imgdur = [];
crossdur = [];
for trial = 1:length(per);
    imgdur = [imgdur per(trial).imgdur];
    crossdur = [crossdur per(trial).crossfixdur];
end
%%
imgs = zeros(2,40);
img_trial = 1;
for trial = 1:length(per);
    cndline = textscan(cndfil(per(trial).cnd+1,:),'%d');
    imgnums = [cndline{1}(end-1) cndline{1}(end)];
    imgs(:,trial) = imgnums';
end

imgnum = zeros(2,40); 
for trial = 1:size(imgs,2)
    itmline1 = textscan(itmfil(imgs(1,trial)+6,:),'%s');
    imgnum1 = str2double(itmline1{1}{end}(11:end-4));
    imgnum(1,trial) = imgnum1;
    
        itmline2 = textscan(itmfil(imgs(2,trial)+6,:),'%s');
    imgnum2 = str2double(itmline2{1}{end}(11:end-4));
    imgnum(2,trial) = imgnum2;
end
%%
LRpos = [];
for trial = 1:size(imgnum,2)
    if imgnum(1,trial) ~= imgnum(2,trial)
        if itmposx(imgnum(2,trial)) < 0
           LRpos = [LRpos -1]; 
        else
            LRpos = [LRpos 1];
        end
    end
end
%%
block_trials = zeros(1,79);
for block = 1:79
    block_trials(block) = length(find(event_arr == 500+block));
end
    
    