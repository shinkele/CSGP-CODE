function [x_algo_cont,x_debias_algo_cont,obj_algo_cont,times_algo_cont,debias_start_algo,mse_algo_cont] = method(algo_name,args)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

% args.y=y;args.R=R;args.tau=tau;args.Debias=debias;args.Continuation=1;
% args.ContinuationSteps=steps;args.FirstTauFactor=first_tau_factor;
% args.Monotone=1;args.Initialization=2;args.True_x=f;args.StopCriterion=stopCri;
% args.ToleranceA=tolA;args.Verbose=0;

 [x_algo_cont,x_debias_algo_cont,obj_algo_cont,times_algo_cont,debias_start_algo,mse_algo_cont]= ...
  feval(algo_name,args.y,args.R,args.tau,...
                  'Debias',args.Debias,...
                  'Continuation',args.Continuation,...
                  'ContinuationSteps',args.ContinuationSteps,...
                  'FirstTauFactor',args.FirstTauFactor,...
                  'Monotone',args.Monotone,...
                  'Initialization',args.Initialization,...
                  'True_x',args.True_x,...
                  'StopCriterion',args.StopCriterion,...
       	          'ToleranceA',args.ToleranceA,...
                  'Verbose',args.Verbose);

end

