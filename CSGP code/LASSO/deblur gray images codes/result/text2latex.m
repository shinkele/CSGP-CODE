clc;clear;

image_name={'baboon.bmp','blobs.bmp','boats.bmp','brain.bmp','brickwall.bmp','bridge.bmp','cameraman.bmp','carpet.bmp',...
    'chart.tiff','circles.tif','clown.bmp','fruits.bmp','girlface.bmp','houses.bmp','kiel.bmp','lena.bmp',...
    'lighthouse.bmp','man.bmp','peppers.bmp','shape.jpg','tank.bmp','tank2.bmp','textureA.bmp',...
    'textureB.bmp','truck.bmp','trucks.bmp','zelda.bmp','zelda2.bmp'};


data=textread('get_result_latex.txt');


fid = fopen('result_latex_version','w');
[row, col]=size(data);

for i=1:row
    temp=data(i,:); % for each row
    
    name=image_name{temp(1)}; % for name
    name(end-3:end)=[];
    
    % iter, time, objective, mse, snr, ssim
%     m1 = temp(i,3:8); 
%     m2 = temp(i,9:14);
%     m3 = temp(i,15:20);
    m1 = temp([3,4,8]); 
    m2 = temp([9,10,14]);
    m3 = temp([15,16,20]);

    str=['%s(%d)'];
    for jj = 1:3
        str=[str,' & %.2f/%.2f/%.3f']; 
    end
    str=[str,' %s \n'];
    fprintf(fid,str,name,temp(2),m1,m2,m3,' \\');
    
end

       
