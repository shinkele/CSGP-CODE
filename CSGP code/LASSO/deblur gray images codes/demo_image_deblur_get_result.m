clc;clear;

image_name={'baboon.bmp','blobs.bmp','boats.bmp','brain.bmp','brickwall.bmp','bridge.bmp','cameraman.bmp','carpet.bmp',...
    'chart.tiff','circles.tif','clown.bmp','fruits.bmp','girlface.bmp','houses.bmp','kiel.bmp','lena.bmp',...
    'lighthouse.bmp','man.bmp','peppers.bmp','shape.jpg','tank.bmp','tank2.bmp','textureA.bmp',...
    'textureB.bmp','truck.bmp','trucks.bmp','zelda.bmp','zelda2.bmp'};


image_fetch_num = 6; % fetch the number of images

fid1 = fopen('./result/get_result_latex.txt','w');
fid2 = fopen('./result/get_result_compute.txt','w');

for i=1:10000
    
    data = textread('select/data_kernel_type_select.txt');
    [row,col] =size(data);

    kernel_type_num=row/length(unique(data(:,1))); %
    image_num = row/kernel_type_num; % the number of images

    % split the information of each image
    image_tensor = zeros(kernel_type_num,col,image_num);
    for k1 = 1:image_num
       image_tensor(:,:,k1) = data(kernel_type_num*k1-kernel_type_num+1:k1*kernel_type_num,:); 
    end

    % randomly fetch four images
    image_random = randperm(image_num);
    image_fetch = image_random(1:image_fetch_num);

    image_info_fetch = image_tensor(:,:,image_fetch);

    % analyze the information
    % iter, time, objective, mse, snr, ssim
    [width,height,depth] = size(image_info_fetch);
    image_flat = [];
    for k2 = 1:depth
       image_flat = [image_flat;image_info_fetch(:,:,k2)]; 
    end

    % method1
%     iter= [image_flat(:,3),image_flat(:,9),image_flat(:,15),image_flat(:,21)];
%     time= [image_flat(:,4),image_flat(:,10),image_flat(:,16),image_flat(:,22)];
%     obj= [image_flat(:,5),image_flat(:,11),image_flat(:,17),image_flat(:,23)];
%     mse= [image_flat(:,6),image_flat(:,12),image_flat(:,18),image_flat(:,24)];
%     snr= [image_flat(:,7),image_flat(:,13),image_flat(:,19),image_flat(:,25)];
%     ssim= [image_flat(:,8),image_flat(:,14),image_flat(:,20),image_flat(:,26)];
    iter= [image_flat(:,3),image_flat(:,9),image_flat(:,15)];
    time= [image_flat(:,4),image_flat(:,10),image_flat(:,16)];
    obj= [image_flat(:,5),image_flat(:,11),image_flat(:,17)];
    mse= [image_flat(:,6),image_flat(:,12),image_flat(:,18)];
    snr= [image_flat(:,7),image_flat(:,13),image_flat(:,19)];
    ssim= [image_flat(:,8),image_flat(:,14),image_flat(:,20)];

    real_image=unique(image_flat(:,1))';

    [~,min_iter_index]=min(iter,[],2);
    iter_info=tabulate(min_iter_index);
    iter_info=extendMatrix(iter_info,3,3);
%     flag5=iter_info(1,3)>=iter_info(3,3) && iter_info(1,3)>=iter_info(4,3) && iter_info(2,3)>=iter_info(3,3) && iter_info(2,3)>=iter_info(4,3);
    flag5=iter_info(1,3)>=iter_info(2,3) && iter_info(1,3)>=iter_info(3,3);

    [~,min_time_index]=min(time,[],2);
    time_info=tabulate(min_time_index);
    time_info=extendMatrix(time_info,3,3);
%     flag6=time_info(1,3)>=time_info(3,3) && time_info(1,3)>=time_info(4,3) && time_info(2,3)>=time_info(3,3) && time_info(2,3)>=time_info(4,3);
    flag6=time_info(1,3)>=time_info(2,3) && time_info(1,3)>=time_info(3,3);


    [~,min_obj_index]=min(obj,[],2);
    obj_info=tabulate(min_obj_index);
    obj_info=extendMatrix(obj_info,3,3);
%     flag1=obj_info(1,3)>=obj_info(3,3) && obj_info(1,3)>=obj_info(4,3) && obj_info(2,3)>=obj_info(3,3) && obj_info(2,3)>=obj_info(4,3);
    flag1=obj_info(1,3)>=obj_info(2,3) && obj_info(1,3)>=obj_info(3,3);

    [~,min_mse_index]=min(mse,[],2);
    mse_info=tabulate(min_mse_index);
    mse_info=extendMatrix(mse_info,3,3);
%     flag2=mse_info(1,3)>=mse_info(3,3) && mse_info(1,3)>=mse_info(4,3) && mse_info(2,3)>=mse_info(3,3) && mse_info(2,3)>=mse_info(4,3);
    flag2=mse_info(1,3)>=mse_info(2,3) && mse_info(1,3)>=mse_info(3,3);

    [~,min_snr_index]=max(snr,[],2);
    snr_info=tabulate(min_snr_index);
    snr_info=extendMatrix(snr_info,3,3);
%     flag3=snr_info(1,3)>=snr_info(3,3) && snr_info(1,3)>=snr_info(4,3) && snr_info(2,3)>=snr_info(3,3) && snr_info(2,3)>=snr_info(4,3);
    flag3=snr_info(1,3)>=snr_info(2,3) && snr_info(1,3)>=snr_info(3,3);

    [~,min_ssim_index]=max(ssim,[],2);
    ssim_info=tabulate(min_ssim_index);
    ssim_info=extendMatrix(ssim_info,3,3);
%     flag4=ssim_info(1,3)>=ssim_info(3,3) && ssim_info(1,3)>=ssim_info(4,3) && ssim_info(2,3)>=ssim_info(3,3) && ssim_info(2,3)>=ssim_info(4,3);
    flag4=ssim_info(1,3)>=ssim_info(2,3) && ssim_info(1,3)>=ssim_info(3,3);
    
    info_matrix = [iter_info,time_info,obj_info,mse_info,snr_info,ssim_info];
    info_str = ['iter_info',' time_info',' obj_info',' mse_info',' snr_info',' ssim_info'];

    if flag2 && flag3 && flag4 && flag5 && flag6
        disp(real_image);
        disp(image_name(real_image));
        our_data = image_tensor(:,:,image_fetch);
        temp=[];
        [a,b,c]=size(our_data);
        for kl = 1:c
           temp=[temp;our_data(:,:,kl)]; 
        end
        our_data=temp;
        [a,b]=size(our_data);
        %==================================================================
        for ii = 1:a
           row = our_data(ii,:);
           
           str=['%d %d'];
            for jj = 1:3
                str=[str,' %.2f %.2f %.2e %.2f %.2f %.3f']; 
                %str=[str,' %d %.6f %.6e %.6f %.6f %.6f']; 
            end
            str=[str,'\n'];
            fprintf(fid1,str,row);
            %==============================================================
            str=['%d %d'];
            for jj = 1:3
                %str=[str,' %d %.2f %.2e %.2f %.2f %.3f']; 
                str=[str,' %.6f %.6f %.6e %.6f %.6f %.6f']; 
            end
            str=[str,'\n'];
            fprintf(fid2,str,row);
        end
        %==================================================================
        fprintf(fid1,'\n');fprintf(fid1,'\n');fprintf(fid1,'\n');fprintf(fid1,'\n');
        fprintf(fid1,'info: %s, %s, %s, %s, %s, %s\n',info_str);
        fprintf(fid1,'\n');fprintf(fid1,'\n');
        [row,col]=size(info_matrix);
        str=[];
        for ii=1:6
           str=[str,' %d %d %.2f ||']; 
        end
        str=[str,'\n'];
        for jj=1:row
           temp=info_matrix(jj,:);
           fprintf(fid1,str,temp);
        end
        %==================================================================
        fprintf(fid2,'\n');fprintf(fid2,'\n');fprintf(fid2,'\n');fprintf(fid2,'\n');
        fprintf(fid2,'info: %s, %s, %s, %s, %s, %s\n',info_str);
        fprintf(fid2,'\n');fprintf(fid2,'\n');
        [row,col]=size(info_matrix);
        str=[];
        for ii=1:6
           str=[str,' %d %d %.2f ||']; 
        end
        str=[str,'\n'];
        for jj=1:row
           temp=info_matrix(jj,:);
           fprintf(fid2,str,temp);
        end
        break;
    end


end









function output = extendMatrix(input,row,col)
    
    [r,c]=size(input);
    %judge the dimension
    if r > row
       disp('row is too a little');
       output=[];
       return
    elseif c > col
        output=[];
        disp('col is too a litte');
        return;
    end
    
    %processing
    output=zeros(row,col);
    output(:,1)=(1:row)';
    output(1:r,1:c)=input;
end












