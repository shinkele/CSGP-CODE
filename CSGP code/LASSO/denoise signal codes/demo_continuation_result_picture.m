% This is a demo for solving \ell_1-norm minimization
% problem
% arg min_x = 0.5*|| y - A x ||_2^2 + tau || x ||_1
% using the algorithm modified three-term conjugate gradient projection method, described in the following paper
%
clc
clear all
close all
addpath Method
randn('seed',1); % 1 for the experiments in the paper
rand('seed',1); % 1 for the experiments in the paper
%%
% parameter
steps = 5;
debias = 0;
stopCri=1;
tolA=1.e-5; % the stopping criterion 
TYPE = 0;
if TYPE == 1
    for index=8 %%%%%%%%%%%%%%%%%%%%%%%%%% 8
        % n is the original signal length
        n=1024*index;
        % k is number of observations to make
        k=256*index;
        % number of spikes to put down
        % n_spikes = floor(.01*n);
        n_spikes = 32*index;  % number of nonzeros
        % generate problems
        %x_algo1_cont=[];

        for repeats=1
            % random +/- 1 signal
            f = zeros(n,1);
            q = randperm(n);
            % f(q(1:n_spikes)) = sign(randn(n_spikes,1));
            f(q(1:n_spikes)) = randn(n_spikes,1);
            disp('Creating measurement matrix...');
            R = randn(k,n);
    %         for j=1:n
    %             R(:,j) = R(:,j)/norm(R(:,j));%normalize R
    %         end
            % orthonormalize rows  
            R = orth(R')'; 
    %         if n==8192
    %             load Rmatrix_2048_8192.mat
    %         else
    %             R = randn(k,n);
                % orthonormalize rows
    %             R = orth(R')'; 
    %         end
    % %         if ~exist('R','var') % measurement matrix
    % %             disp('Creating measurement matrix...');
    % %             R = randn(k,n);
    % %             % orthonormalize rows
    % %             R = orth(R')';
    % %         end
    % %         if n == 8192  % in this case, we load a precomputed matrix to save some time
    % %             load Rmatrix_2048_8192.mat
    % %         end
            disp('Finished creating matrix');   
            %% Define all its functions:
            % Create handles to functions R(x)=Rx and R^{T}(x)=R^{T}x:
            hR = @(x) R*x;
            hRt = @(x) (x'*R)';% R'*x;
            % noisy observations
            sigma = 0.0001;  % sigma = 0.0001;
            y = hR(f) + sigma*randn(k,1);
            % regularization parameter
            tau = 0.005*max(abs(R'*y));
            first_tau_factor = 0.8*(max(abs(R'*y))/tau);
            %%
            %==================================================================
            args.y=y;args.R=R;args.tau=tau;args.Debias=debias;args.Continuation=1;
            args.ContinuationSteps=steps;args.FirstTauFactor=first_tau_factor;
            args.Monotone=1;args.Initialization=2;args.True_x=f;args.StopCriterion=stopCri;
            args.ToleranceA=tolA;args.Verbose=0;
            %==================================================================
            method_name={'LCGP1_CS','HSDY_CS','CGD_CS'};
            %%
            %==================================================================

            disp(['Starting',' ',method_name{1},' ','monotonic with continuation'])
            [x_algo1_cont,x_debias_algo1_cont,obj_algo1_cont,times_algo1_cont,debias_start_algo1,mse_algo1_cont]=method(method_name{1},args);
            t_algo1_cont = times_algo1_cont(end);
            Itr_algo1=length(obj_algo1_cont);
            Obj_algo1 = obj_algo1_cont(end);
            MSE_algo1 = (1/n)*norm(x_algo1_cont-f)^2;

            disp(['Starting',' ',method_name{2},' ','monotonic with continuation'])
            [x_algo2_cont,x_debias_algo2_cont,obj_algo2_cont,times_algo2_cont,debias_start_algo2,mse_algo2_cont]=method(method_name{2},args);
            t_algo2_cont = times_algo2_cont(end);
            Itr_algo2=length(obj_algo2_cont);
            Obj_algo2 = obj_algo2_cont(end);
            MSE_algo2 = (1/n)*norm(x_algo2_cont-f)^2;

            disp(['Starting',' ',method_name{3},' ','monotonic with continuation'])
            [x_algo3_cont,x_debias_algo3_cont,obj_algo3_cont,times_algo3_cont,debias_start_algo3,mse_algo3_cont]=method(method_name{3},args);
            t_algo3_cont = times_algo3_cont(end);
            Itr_algo3=length(obj_algo3_cont);
            Obj_algo3 = obj_algo3_cont(end);
            MSE_algo3 = (1/n)*norm(x_algo3_cont-f)^2;

%             disp(['Starting',' ',method_name{4},' ','monotonic with continuation'])
%             [x_algo4_cont,x_debias_algo4_cont,obj_algo4_cont,times_algo4_cont,debias_start_algo4,mse_algo4_cont]=method(method_name{4},args);
%             t_algo4_cont = times_algo4_cont(end);
%             Itr_algo4=length(obj_algo4_cont);
%             Obj_algo4 = obj_algo4_cont(end);
%             MSE_algo4 = (1/n)*norm(x_algo4_cont-f)^2;

            % =============================================================
            save('./result/picture_data/n.mat','n');
            save('./result/picture_data/k.mat','k');
            save('./result/picture_data/tau.mat','tau');

            save('./result/picture_data/mse_algo1_cont.mat','mse_algo1_cont');
            save('./result/picture_data/mse_algo2_cont.mat','mse_algo2_cont');
            save('./result/picture_data/mse_algo3_cont.mat','mse_algo3_cont');
%             save('./result/picture_data/mse_algo4_cont.mat','mse_algo4_cont');

            save('./result/picture_data/times_algo1_cont.mat','times_algo1_cont');
            save('./result/picture_data/times_algo2_cont.mat','times_algo2_cont');
            save('./result/picture_data/times_algo3_cont.mat','times_algo3_cont');
%             save('./result/picture_data/times_algo4_cont.mat','times_algo4_cont');
            
            save('./result/picture_data/obj_algo1_cont.mat','obj_algo1_cont');
            save('./result/picture_data/obj_algo2_cont.mat','obj_algo2_cont');
            save('./result/picture_data/obj_algo3_cont.mat','obj_algo3_cont');
%             save('./result/picture_data/obj_algo4_cont.mat','obj_algo4_cont');

            save('./result/picture_data/times_algo1_cont.mat','times_algo1_cont');
            save('./result/picture_data/times_algo2_cont.mat','times_algo2_cont');
            save('./result/picture_data/times_algo3_cont.mat','times_algo3_cont');
%             save('./result/picture_data/times_algo4_cont.mat','times_algo4_cont');
            
            save('./result/picture_data/f.mat','f');
            save('./result/picture_data/n_spikes.mat','n_spikes');
            save('./result/picture_data/y.mat','y');
            save('./result/picture_data/x_algo1_cont.mat','x_algo1_cont');
            save('./result/picture_data/x_algo2_cont.mat','x_algo2_cont');
            save('./result/picture_data/x_algo3_cont.mat','x_algo3_cont');
%             save('./result/picture_data/x_algo4_cont.mat','x_algo4_cont');
            
            
        end
    end
end
%==========================================================================



%================= Plotting results ==========obj_IST,times_IST,debias_s,mses_IST
method_name={'CSGP','HSDY','CGD'};
% method_name{1}(end-2:end)=[];
% method_name{2}(end-2:end)=[];
% method_name{3}(end-2:end)=[];
%method_name={'CSGP','HSDY','CGD'};
%currentFolder = pwd;

load('./result/picture_data/n.mat');
load('./result/picture_data/k.mat');
load('./result/picture_data/tau.mat');

load('./result/picture_data/mse_algo1_cont.mat');
load('./result/picture_data/mse_algo2_cont.mat');
load('./result/picture_data/mse_algo3_cont.mat');
% load('./result/picture_data/mse_algo4_cont.mat');

load('./result/picture_data/times_algo1_cont.mat');
load('./result/picture_data/times_algo2_cont.mat');
load('./result/picture_data/times_algo3_cont.mat');
% load('./result/picture_data/times_algo4_cont.mat');

load('./result/picture_data/obj_algo1_cont.mat');
load('./result/picture_data/obj_algo2_cont.mat');
load('./result/picture_data/obj_algo3_cont.mat');
% load('./result/picture_data/obj_algo4_cont.mat');

load('./result/picture_data/times_algo1_cont.mat');
load('./result/picture_data/times_algo2_cont.mat');
load('./result/picture_data/times_algo3_cont.mat');
% load('./result/picture_data/times_algo4_cont.mat');

load('./result/picture_data/f.mat');
load('./result/picture_data/n_spikes.mat');
load('./result/picture_data/y.mat');
load('./result/picture_data/x_algo1_cont.mat');
load('./result/picture_data/x_algo2_cont.mat');
load('./result/picture_data/x_algo3_cont.mat');
% load('./result/picture_data/x_algo4_cont.mat');

figure(1)
plot(mse_algo1_cont,'b:','LineWidth',2)
hold on
plot(mse_algo2_cont,'r-','LineWidth',2)
hold on
plot(mse_algo3_cont,'g-','LineWidth',2)
hold on
% plot(mse_algo4_cont,'k-','LineWidth',2)
% hold on
legend(method_name);
set(gca,'FontName','Times','FontSize',16)
xlabel('Iterations')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
axis on
grid on
hold off

figure(2)
plot(times_algo1_cont,mse_algo1_cont,'b:','LineWidth',2)
hold on
plot(times_algo2_cont,mse_algo2_cont,'r-','LineWidth',2)
hold on
plot(times_algo3_cont,mse_algo3_cont,'g-','LineWidth',2)
hold on
% plot(times_algo4_cont,mse_algo4_cont,'k-','LineWidth',2)
% hold on
legend(method_name)
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
axis on
grid on
hold off
%==========================================================================


figure(3)
plot(obj_algo1_cont,'b:','LineWidth',2)
hold on
plot(obj_algo2_cont,'r-','LineWidth',2)
hold on
plot(obj_algo3_cont,'g-','LineWidth',2)
hold on
% plot(obj_algo4_cont,'k-','LineWidth',2)
legend(method_name)
set(gca,'FontName','Times','FontSize',16)
xlabel('Iterations')
ylabel('ObjFun')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
axis on
grid on
hold off

figure(4)
plot(times_algo1_cont,obj_algo1_cont,'b:','LineWidth',2)
hold on
plot(times_algo2_cont,obj_algo2_cont,'r-','LineWidth',2)
hold on
hold on
plot(times_algo3_cont,obj_algo3_cont,'g-','LineWidth',2)
hold on
% plot(times_algo4_cont,obj_algo4_cont,'k-','LineWidth',2)
% hold on
legend(method_name)
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('ObjFun')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
axis on
grid on
hold off
%==========================================================================



figure(5)
%scrsz = get(0,'ScreenSize');
%set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(5,1,1)
plot(f,'LineWidth',1.1)
top = max(f(:));
bottom = min(f(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Original (n = %g, number of nonzeros = %g)',n,n_spikes))
axis(v)

% scrsz = get(0,'ScreenSize');
% set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(5,1,2)
plot(y,'LineWidth',1.1)
top = max(y(:));
bottom = min(y(:));
v = [0 k+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Measurement'))
axis(v)

%scrsz = get(0,'ScreenSize');
%set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(5,1,3)
plot(x_algo1_cont(:),'LineWidth',1.1)
top = max(x_algo1_cont(:));
bottom = min(x_algo1_cont(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf([method_name{1},'(MSE = %5.2e, Iter=%g, Time=%4.2fs)'],  (1/n)*norm(x_algo1_cont-f)^2,length(times_algo1_cont),times_algo1_cont(end)))
axis(v)

%scrsz = get(0,'ScreenSize');
%set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(5,1,4)
plot(x_algo2_cont(:),'LineWidth',1.1)
top = max(x_algo2_cont(:));
bottom = min(x_algo2_cont(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf([method_name{2},'(MSE = %5.2e,Iter=%g, Time=%4.2fs)'],(1/n)*norm(x_algo2_cont-f)^2,length(times_algo2_cont),times_algo2_cont(end)))
axis(v)

%scrsz = get(0,'ScreenSize');
%set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(5,1,5)
plot(x_algo3_cont(:),'LineWidth',1.1)
top = max(x_algo3_cont(:));
bottom = min(x_algo3_cont(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf([method_name{3},'(MSE = %5.2e,Iter=%g, Time=%4.2fs)'],(1/n)*norm(x_algo3_cont-f)^2,length(times_algo3_cont),times_algo3_cont(end)))
axis(v)

%scrsz = get(0,'ScreenSize');
%set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
% subplot(6,1,6)
% plot(x_algo4_cont(:),'LineWidth',1.1)
% top = max(x_algo4_cont(:));
% bottom = min(x_algo4_cont(:));
% v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
% set(gca,'FontName','Times')
% set(gca,'FontSize',14)
% title(sprintf([method_name{4},'(MSE = %5.2e,Iter=%g, Time=%4.2fs)'],(1/n)*norm(x_algo4_cont-f)^2,length(times_algo4_cont),times_algo4_cont(end)))
% axis(v)



