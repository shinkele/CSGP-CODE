clear all
clc;
format long
format compact
%% 
% set parameter
randn('seed',1); % 1 for the experiments in the paper
rand('seed',1); % 1 for the experiments in the paper
addpath Method
%%
N =[1 2 3 4 5 6 7 8];
% the dimension of varaibles
n = [1000 5000 10000 50000 100000];
% the set of initial point
I = 1:7;
% the set of algorithm
algorithm_set = {'LCGP1_CS','HSDY_CS','CGD_CS'};
%%
TYPE = 1;
if TYPE == 0
    fid = fopen('result/data.txt','w');  
    problem_set = N;
    for i = problem_set % problems
        for j = 1:length(n) % dimension
            dim = n(j); 
            for k = I 
                repeat_table = [];
                for repeat=1:10 % to set a number
                    x0 = rand(dim,1); % generate a random initial points 
                    row_table = [];
                    for l = 1:length(algorithm_set) % algorithm
                        fprintf('Problem: %d, Dim: %d, repeat: %d, algorithm: %s\n',i,dim,k,algorithm_set{l});
                        [Tcpu,NF,Itr,NG] = feval(algorithm_set{l},dim,x0,i); 
                        row_table = [row_table,Tcpu,NF,Itr,NG];
                    end
                    repeat_table =[repeat_table;row_table];
                end
                mean_table=mean(repeat_table);
                %fprintf(fid,'%d %d %d %.2e %.2f %.2f %.2e %.2e %.2f %.2f %.2e %.2e %.2f %.2f %.2e\n',i,dim,k,mean_table);
                str = ['%d %d %d'];
                for kk = 1:length(algorithm_set)
                    str = [str,' ','%.2e %2.f %2.f %.2e']; 
                end
                str = [str,'\n'];
                fprintf(fid,str,i,dim,k,mean_table);
            end
        end
    end
    fclose(fid);
end
%%
data=textread('result/data.txt'); % read the text file
[row,col]=size(data);

TIME = data(:,4:4:col);
ITER = data(:,6:4:col);
NFF = data(:,5:4:col);

figure(1);
perf(TIME,'logplot');
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend(algorithm_set,'Location','SouthEast');   
saveas(gcf,['result/CPU.jpg']);
hold off

figure(2);
perf(NFF,'logplot');
xlabel('\tau','Interpreter','tex');   
ylabel('\rho(\tau)','Interpreter','tex');               
legend(algorithm_set,'Location','SouthEast'); 
saveas(gcf,['result/NFF.jpg']);
hold off

figure(3);
perf(ITER,'logplot');
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend(algorithm_set,'Location','SouthEast'); 
saveas(gcf,['result/NI.jpg']);
hold off





