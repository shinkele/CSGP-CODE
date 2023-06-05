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
fid=fopen('result/data_result.txt','w');
fprintf(fid,'%s & %s & %s & %s\n','k','n','n_spikes','algo(iter,time,MSE)');
for index=1:0.5:10
    % n is the original signal length
    n=1024*index;
    % k is number of observations to make
    k=256*index;
    % number of spikes to put down
    % n_spikes = floor(.01*n);
    n_spikes = 32*index;  % number of nonzeros
    % generate problems
    progress_r=[];
    for repeats=1:10
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
        method_name={'LCGP1_CS','HSDY_CS','CGD_CS'};
        %%
        %==================================================================
        message=[];
        for ii = 1:length(method_name)
            disp(['Starting',' ',method_name{ii},' ','monotonic with continuation'])
            [x_algo_cont,x_debias_algo_cont,obj_algo_cont,times_algo_cont,debias_start_algo,mse_algo_cont]=method(method_name{ii},args);
            t_algo_cont = times_algo_cont(end);
            Itr_algo=length(obj_algo_cont);
            %Obj_CGDCS = obj_CGDCS_cont(end);
            MSE_algo = (1/n)*norm(x_algo_cont-f)^2;
            message=[message Itr_algo t_algo_cont MSE_algo];
        end

        progress_r=[progress_r;message];
    end
    fprintf(fid,'(%d, %d, %d) & %.2f/%.2f/%.3e & %.2f/%.2f/%.3e & %.2f/%.2f/%.3e%s\n', ... 
        k,n,n_spikes,mean(progress_r,1),' \\');
end
fclose(fid);

% %
% fprintf(1,'\n\n-------------------------------------------------\n')   
% fprintf(1,'-------------------------------------------------\n')   
% fprintf(1,'Problem: n = %g,  k = %g, number of spikes = %g\n',n,k,n_spikes)
% fprintf(1,'Parameters: sigma = %g, tau = %g, debiasing = %g\n',sigma,tau,debias)
% fprintf(1,'All algorithms initialized with zeros\n')
% fprintf(1,'-------------------------------------------------\n')
% 
% fprintf(1,'\nCGDESCENT_CS monotone continuation; cpu: %6.2f secs (%d iterations)\n',...
%         t_CGDCS_cont,length(obj_CGDCS_cont))
% fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
%           obj_CGDCS_cont(end),(1/n)*norm(x_CGDCS_cont-f)^2)      
% 
% fprintf(1,'\nHTTCGP_CS monotone continuation; cpu: %6.2f secs (%d iterations)\n',...
%         t_SGCS_cont,length(obj_SGCS_cont))
% fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
%           obj_SGCS_cont(end),(1/n)*norm(x_SGCS_cont-f)^2)      
% 
%       
% fprintf(1,'-------------------------------------------------\n')
% fprintf(1,'-------------------------------------------------\n')
% 
% % ================= Plotting results ==========obj_IST,times_IST,debias_s,mses_IST
% %%
% figure(1)
% plot(mse_SGCS_cont,'b:','LineWidth',2)
% hold on
% plot(mse_CGDCS_cont,'r-','LineWidth',2)
% hold on
% legend('HTTCGP-CS','CGD') % SG-CS
% set(gca,'FontName','Times','FontSize',16)
% xlabel('Iterations')
% ylabel('MSE')
% title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
% axis on
% grid on
% hold off
% 
% figure(2)
% plot(times_SGCS_cont,mse_SGCS_cont,'b:','LineWidth',2)
% hold on
% plot(times_CGDCS_cont,mse_CGDCS_cont,'r-','LineWidth',2)
% legend('HTTCGP-CS','CGD')
% set(gca,'FontName','Times','FontSize',16)
% xlabel('CPU time (seconds)')
% ylabel('MSE')
% title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
% axis on
% grid on
% hold off
% % %==========================================================================
% figure(3)
% plot(obj_SGCS_cont,'b:','LineWidth',2)
% hold on
% plot(obj_CGDCS_cont,'r-','LineWidth',2)
% legend('HTTCGP-CS','CGD')
% set(gca,'FontName','Times','FontSize',16)
% xlabel('Iterations')
% ylabel('ObjFun')
% title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
% axis on
% grid on
% hold off
% % 
% figure(4)
% plot(times_SGCS_cont,obj_SGCS_cont,'b:','LineWidth',2)
% hold on
% plot(times_CGDCS_cont,obj_CGDCS_cont,'r-','LineWidth',2)
% hold on
% legend('HTTCGP-CS','CGD')
% set(gca,'FontName','Times','FontSize',16)
% xlabel('CPU time (seconds)')
% ylabel('ObjFun')
% title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
% axis on
% grid on
% hold off
% %==========================================================================
% %%
% figure(5)
% scrsz = get(0,'ScreenSize');
% % set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
% subplot(4,1,1)
% plot(f,'LineWidth',1.1)
% top = max(f(:));
% bottom = min(f(:));
% v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
% set(gca,'FontName','Times')
% set(gca,'FontSize',14)
% title(sprintf('Original (n = %g, number of nonzeros = %g)',n,n_spikes))
% axis(v)
% 
% scrsz = get(0,'ScreenSize');
% % set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
% subplot(4,1,2)
% plot(y,'LineWidth',1.1)
% top = max(y(:));
% bottom = min(y(:));
% v = [0 k+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
% set(gca,'FontName','Times')
% set(gca,'FontSize',14)
% title(sprintf('Measurement'))
% axis(v)
% 
% scrsz = get(0,'ScreenSize');
% % set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
% subplot(4,1,3)
% plot(x_SGCS_cont(:),'LineWidth',1.1)
% top = max(x_SGCS_cont(:));
% bottom = min(x_SGCS_cont(:));
% v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
% set(gca,'FontName','Times')
% set(gca,'FontSize',14)
% title(sprintf('HTTCGP-CS (MSE = %5.2e, Iter=%g, Time=%4.2fs)',  (1/n)*norm(x_SGCS_cont-f)^2,length(times_SGCS_cont),times_SGCS_cont(end)))
% axis(v)
% 
% scrsz = get(0,'ScreenSize');
% % set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
% subplot(4,1,4)
% plot(x_CGDCS_cont(:),'LineWidth',1.1)
% top = max(x_CGDCS_cont(:));
% bottom = min(x_CGDCS_cont(:));
% v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
% set(gca,'FontName','Times')
% set(gca,'FontSize',14)
% title(sprintf('CGD (MSE = %5.2e,Iter=%g, Time=%4.2fs)',(1/n)*norm(x_CGDCS_cont-f)^2,length(times_CGDCS_cont),times_CGDCS_cont(end)))
% axis(v)
% % 
% % 
