clc;clear;


for kernel_type = 1:4

load(['./picture_data/recimg_kernel_',num2str(kernel_type),'.mat'],'recimg');
load(['./picture_data/blurred_kernel_',num2str(kernel_type),'.mat'],'blurred');
load(['./picture_data/original_kernel_',num2str(kernel_type),'.mat'],'original');
[row,col]=size(recimg);

for i=1:row

    figure(1)
    imagesc(original{i})
    colormap(gray(255))
    axis off
    axis equal
    title('Original image','FontName','Times','FontSize',22)

    figure(2)
    imagesc(blurred{i})
    colormap(gray(255))
    axis off
    axis equal
    title('Blurred image','FontName','Times','FontSize',22)

    for j=1:col
        methods_name={'CSGP','HSDY','CGD'};
        figure(j+2)
        imagesc(recimg{i,j})
        colormap(gray(255))
        axis off
        axis equal
        title(methods_name{j},'FontName','Times','FontSize',22)
    end
    pause();
end
    


end 






















