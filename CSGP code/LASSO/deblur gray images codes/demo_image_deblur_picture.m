% This demo shows the reconstruction of a sparse image
% of randomly placed plus an minus ones
% 
%%
clc
clear all
close all
%clf
randn('seed',1); % 1 for the experiments in the paper
rand('seed',1); % 1 for the experiments in the paper
addpath Images
addpath Method
addpath Core
%%
% the test images set
image_name={'baboon.bmp','blobs.bmp','boats.bmp','brain.bmp','brickwall.bmp','bridge.bmp','cameraman.bmp','carpet.bmp',...
    'chart.tiff','circles.tif','clown.bmp','fruits.bmp','girlface.bmp','houses.bmp','kiel.bmp','lena.bmp',...
    'lighthouse.bmp','man.bmp','peppers.bmp','shape.jpg','tank.bmp','tank2.bmp','textureA.bmp',...
    'textureB.bmp','truck.bmp','trucks.bmp','zelda.bmp','zelda2.bmp'};


TYPE = 1;
verbose_output = 1;
verbose_image = 1; % 1: plot the picture
repeat_num=1;
for kernel_type=1:4 % choose the type of kernel: 1,2,3,4
image_select =13; 
Gausian_std = 2.25;  
%%
if TYPE == 1
    recimg=cell(1,3);
    original=cell(1,1);
    blurred=cell(1,1);
    % iter, time, objective, mse, snr, ssim
    iter_total=[];time_total=[];obj_total=[];mse_total=[];snr_total=[];ssim_total=[];
    %fid_data=fopen('result/data_Gausian_kernel_picture.txt','w');
    %fid_image=fopen('select/image_Gausian_kernel_select.txt','w');
    pic_counter = 0;
    for pic = image_select
        im=image_name{pic};
        pic_counter = pic_counter + 1;

        f=double(imread(im)); 
        [width,height] = size(f);

        % to test whether the width and height of image is even.
        if mod(width,2) ==  1
            f(1,:,:)=[];
        end
        if mod(height,2) == 1
            f(:,1,:)=[];
        end
        [width,height] = size(f);

        if verbose_image == 1
            figure(1)
            imagesc(f)
            original{pic_counter,1}=f;
            colormap(gray(255))
            axis off
            axis equal
            title('Original image','FontName','Times','FontSize',22)
        end
        %%
        term_noise = randn(size(f)); % to add the noise
        gausian_message=[]; 
        for g1 = Gausian_std
            % create observation operator; in this case 
            % it will be a blur function composed with an
            % inverse weavelet transform
            disp('Creating observation operator...');

            middle2 = height/2 + 1;
            middle1 = width/2 + 1;

            switch kernel_type
                case 1
                    %uncomment the following lines for Experiment 1 (see paper)
                    sigma = sqrt(0.56);
                    h = zeros(size(f));
                    for i=-4:4
                       for j=-4:4
                          h(i+middle1,j+middle2)= 1; 
                       end
                    end
                case 2
                    % uncomment the following lines for Experiment 2 (see paper)
                    sigma = sqrt(2);
                    %sigma = sqrt(0.01);
                    h = zeros(size(f)); % h is the same size as f and all zeros.
                    for i=-4:4
                       for j=-4:4
                          h(i+middle1,j+middle2)= (1/(1+i*i+j*j));
                       end
                    end
                case 3
                    % uncomment the following lines for Experiment 2 (see paper)
                    sigma = sqrt(4);
                    %sigma = sqrt(0.01);
                    h = zeros(size(f)); % h is the same size as f and all zeros.
                    for i=-4:4
                       for j=-4:4
                          h(i+middle1,j+middle2)= (1/(1+i*i+j*j));
                       end
                    end
                case 4
                    % Gaussian kernel
                    sigma = sqrt(8);
                    h = zeros(size(f));
                    for i=-4:4
                       for j=-4:4
                          h(i+middle1,j+middle2)= (1/(2*pi*g1^2))*exp(-(i^2+j^2)/(2*g1^2));
                       end
                    end
            end % for switch

            % % center and normalize the blur
            h = fftshift(h);  % ��Ƶ��ͼ��λ����Ƶ����Ƶ��ͼ����
            % fftshift is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
            h = h/sum(h(:));  % sum(h(:)): sum all elements of h.

            % definde the function handles that compute 
            % the blur and the conjugate blur.
            R = @(x) real(ifft2(fft2(h).*fft2(x))); % fft2(X) returns the two-dimensional Fourier transform of matrix X.
            % ifft2(F) returns the two-dimensional inverse Fourier transform of matrix F
            RT = @(x) real(ifft2(conj(fft2(h)).*fft2(x)));

            % define the function handles that compute 
            % the products by W (inverse DWT) and W' (DWT)
            wav = daubcqf(2);
            for L=3:-1:1
                if mod(width,2^L) == 0 && mod(height,2^L) == 0
                   break; 
                end
            end
            W = @(x) midwt(x,wav,3);
            WT = @(x) mdwt(x,wav,3);

            % Finally define the function handles that compute 
            % the products by A = RW  and A' =W'*R' 
            A = @(x) R(W(x));
            AT = @(x) WT(RT(x));

            repeats_message=[];
            for repeat = 1:repeat_num
                % generate noisy blurred observations
                y = R(f) + sigma*term_noise;

                if verbose_image == 1
                    figure(2)
                    imagesc(y)
                    blurred{pic_counter,1}=y;
                    colormap(gray(255))
                    axis off
                    axis equal
                    title('Blurred image','FontName','Times','FontSize',22)
                end

                %%
                % regularization parameter
                tau = .35;

                % set tolA
                tolA = 1.e-5;

                %==========================================================================
                args.y=y;args.A=A;args.tau=tau;args.Debias=0;args.AT=AT;
                args.True_x=WT(f);args.Initialization=AT(y);
                args.StopCriterion=1;args.ToleranceA=tolA;args.Verbose=0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                methods_name={'LCGP1_CS','HSDY_CS','CGD_CS'};
                
                method_message=[];
                for k1 = 1:length(methods_name)
                    disp(['Starting algorithm',' ',methods_name{k1},' ','to recover the image ...']);

                    [x_algo,theta_debias,obj_algo,t_algo,debias_algo,mses_algo] = method(methods_name{k1},args);

                    SNR_algo=20*log10(norm(f,'fro')/norm(f-W(x_algo),'fro'));
                    ssim_value = ssim(f,W(x_algo));
                    % iter, time, objective, mse, snr, ssim
                    method_message=[method_message,length(t_algo),t_algo(end),obj_algo(end),mses_algo(end)/width/height,SNR_algo,ssim_value];

                    if verbose_output == 1
                        fprintf('%s, %.2f, %s, iter=%d, time=%.2f, obj=%.2e, MSE=%.2f, SNR=%.2f, SSIM=%.3f\n',im,g1,methods_name{k1},length(t_algo),...
                        t_algo(end),obj_algo(end),mses_algo(end)/width/height,SNR_algo,ssim_value);
                    end
                    if verbose_image == 1
                        figure(k1+2)
                        
                        recimg{pic_counter,k1}=W(x_algo);
                        
                        imagesc(W(x_algo))
                        colormap(gray(255))
                        axis off
                        axis equal
                        title(methods_name{k1},'FontName','Times','FontSize',22)
                    end
                    
                end % for method   
                                
                repeats_message=[repeats_message;method_message];
                repeat_message = mean(repeats_message,1);
            end % for repeat
%             str=['%d %.2f'];
%             for tt = 1:length(methods_name)
%                str=[str,' %d %.2f %.2e %.2f %.2f %.3f']; 
%             end
%             str=[str,'\n'];
%             fprintf(fid_data,str,pic,g1,repeat_message);
            gausian_message=[gausian_message;repeat_message];
        end % for guasian
    end % for image
    %fclose(fid_data);
    save(['./result/picture_data/recimg_kernel_',num2str(kernel_type),'.mat'],'recimg');
    save(['./result/picture_data/blurred_kernel_',num2str(kernel_type),'.mat'],'blurred');
    save(['./result/picture_data/original_kernel_',num2str(kernel_type),'.mat'],'original');
% else
% 
%     load('./result/recimg.mat');
%     load('./result/blurred.mat');
%     load('./result/original.mat');
%     [row,col]=size(recimg);
%     
%     for i=1:row
%         
%         figure(1)
%         imagesc(original{i})
%         colormap(gray(255))
%         axis off
%         axis equal
%         title('Original image','FontName','Times','FontSize',22)
% 
%         figure(2)
%         imagesc(blurred{i})
%         colormap(gray(255))
%         axis off
%         axis equal
%         title('Blurred image','FontName','Times','FontSize',22)
%         
%         for j=1:col;
%             figure(j+2)
%             imagesc(recimg{i,j})
%             colormap(gray(255))
%             axis off
%             axis equal
%             title(methods_name{j},'FontName','Times','FontSize',22)
%         end
% %        pause();
%     end
    
end


end

























