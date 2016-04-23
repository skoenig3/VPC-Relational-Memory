main_dir = 'R:\Buffalo Lab\eblab\Cortex Programs\VPC Relational Memory Project\';
sub_dir = {'2006 (Yr1)VPC Images\','2007 (Yr2)VPC Images\','2008 (Yr3)VPC Images\'};
img_prefix = {'I ','III ','IV '};


for s = 1:length(sub_dir)
    
    new_dir = [main_dir sub_dir{s}(1:end-1) '_Cut\'];
    mkdir(new_dir)
    
    for img = 1:20
        if img < 10
            nov_name = [img_prefix{s} '0' num2str(img)];
        else
            nov_name = [img_prefix{s} num2str(img)];
        end
        
        nov_img = imread([main_dir sub_dir{s} nov_name '.bmp']);
        
        nov_img = nov_img(6:end-6,7:259,:);
        nov_img = imresize(nov_img,[192 192]);
        
        imwrite(nov_img,[new_dir num2str(2*img-1) '.bmp'],'bmp')
        
        if exist([main_dir sub_dir{s} nov_name 'L.bmp'])
            LR = 'L';
            rep_img = imread([main_dir sub_dir{s} nov_name 'L.bmp']);
        elseif exist([main_dir sub_dir{s} nov_name 'R.bmp'])
            LR = 'R';
            rep_img = imread([main_dir sub_dir{s} nov_name 'R.bmp']);
        else
            error('Cant find second image')
        end
        
        if LR == 'L'
            rep_img = rep_img(6:end-6,7:259,:);
        else
            rep_img = rep_img(6:end-6,484:end-6,:);
        end
        
          rep_img = imresize(rep_img,[192 192]);
        if LR == 'L'
            imwrite(rep_img,[new_dir num2str(2*img) 'L.bmp'],'bmp')
        else
            imwrite(rep_img,[new_dir num2str(2*img) 'R.bmp'],'bmp')
        end
    end
end

