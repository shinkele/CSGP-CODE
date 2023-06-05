
 function [Tcpu,NF,Itr,NG] = CGD_CS(n,x0,nprob) 
format long
Step=0:4;  % ����һ��5ά��������0,1,2,3,4��
nextstep=Step(1);
% disp(date);
finish=0;
Itr=0;  %% ����
NF=1;   %% Ŀ�꺯���������
%NP=1;  %% ͶӰ�������
tic     %% ��ʼ��ʱ��
while finish==0
    loop=1;
    while loop==1
        switch nextstep
            case Step(1)
                % k=1;
                %% Step 0 ��ʼ�� ��������
                %% ���ó�ʼ��
                epsilon=1e-6;
                epsilon1=1e-7;
                %% ����������
                sigma=0.0001; % sigma=0.01
                gamma=1;    % gamma=1, initial guess for the steplength
                rho=0.1;    % rho=0.9
                %% ��ʼ����
                Fk=problem(n,x0,nprob);  
                dk=(-1)*Fk;
                nextstep=Step(2);
%               tic;
            case Step(2)
                %% Step 1 ��ֹ׼���������
                if norm(Fk,2)<=epsilon || Itr>2000  %����ѭ��
                    %                     disp('��ϲ���������ﵽ����Ҫ��')
                    %                     disp('xk��KKT�㣡��')
                    nextstep=Step(5);  % ����
                    break;
                else  % ������
                 %% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    t =gamma;  % initial steplength
                    L1=1;
                    while L1==1
                        z_new=x0+t*dk;
                        Fz_new=problem(n,z_new,nprob);
                        NF=NF+1;
                        % eta_k=max(0.001,min(0.8,norm(Fz_new,2)));  % 0.001,8
                      %% check the Armijo-type line search condition
                        if (-Fz_new'*dk < sigma*t*norm(dk,2)^2 && t>10^(-10))  % the Armijo-type line search condition violated
                            L1=1;
                            t=t*rho;
                        else
                            L1=0;
                        end
                    end       % ��ֹwhile L1==1
                 %% End Armijo-type line search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     t
                    zk=x0+t*dk;  % zk=z_new;
                    Fzk= problem(n,zk,nprob);  %  Fzk=Fz_new;
%                     if norm(Fzk,2)<=epsilon
%                        nextstep=Step(5); % ����
%                        break;
%                     else
%                        nextstep=Step(3);
%                     end
                end  % ��ֹ if norm(Fk,2)<=epsilon || Itr>2000  
                nextstep=Step(3);
            case Step(3) % ͶӰ
                xik=Fzk'*(x0-zk)/norm(Fzk)^2;
                zk1=x0 - 1.6 * xik * Fzk;     % �˴�ͨ����ͶӰ��ϵ����1
                x1=projectionTechnique(n,zk1,nprob);     % ͶӰ�õ��µĵ�����x_{k+1}
                %NP=NP+1;
                Fk0=Fk;
                Fk=problem(n,x1,nprob);   % �����µĺ���ֵF_{k+1}
                NF=NF+1;
                nextstep=Step(4);
            case Step(4)
                %% Step 2 ����
                Itr=Itr+1;
%               Fk0=feval(nprob,n,x0,1);
%                 sk=x1-x0;
%                 yk=Fk-Fk0;
                %==========================================================
%                 switch method
%                     case 'WYL'
%                         betak=Fk'*(Fk-norm(Fk)/norm(Fk0)*Fk0)/norm(Fk0)^2; % WYL
%                         dk=-1*Fk+betak*dk; 
%                     case 'PRP+'
%                         betak=max(0,Fk'*(Fk-Fk0)/norm(Fk0)^2);
%                         dk=-1*Fk+betak*dk; 
%                    case'HTTCG'
%                         tk=min(0.3,max(0,1-yk'*sk/norm(yk)^2)); % yk1=yk+0.01*dk����ֵЧ��û��ֱ��ʹ��ykЧ����
%                         mu=0.2;   
%                         fenmu=max(norm(Fk0)^2,max(mu*norm(dk)*norm(yk),dk'*yk));
%                         thetak=tk*Fk'*dk/fenmu;  % ԭʼ�ķ�ĸ�� norm(Fk0)^2  
%                         betak=Fk'*yk/fenmu-norm(yk)^2*Fk'*dk/fenmu^2; %(max(mu*norm(dk)*norm(Fk),norm(Fk0)^2))
%                         dk=-Fk+betak*dk+thetak*yk;
%                 end                  
                  gammak = Fk - Fk0;
                  lambdak = 1+norm(Fk)^(-1)*max(0,-gammak'*dk*t/norm(t*dk)^2);
                  yk=gammak+lambdak*t*norm(Fk)*dk;
                  betak=(yk'*Fk-2*dk'*Fk*norm(yk)^2/(dk'*yk))/(dk'*yk);
                  dk=-Fk+betak*dk;
                %==========================================================
                x0=x1; 
                if norm(dk)<epsilon1
                    nextstep=Step(5);
                    break;
                end
%                 t=t*(Fk11'*dk0)/(Fk'*dk);
%                 if (Fk'*dk>epsilon1)
%                     x1=x0+t*(-Fk);
%                 end
                nextstep=Step(2);
            case Step(5)
                loop=0;
%                 finish=1;
%                 break;
        end
    end
    if Itr>=2001 %|| norm(gk,2)>epsilon || isnan(Itr)==1  || abs(t)<epsilon  norm(gk,2)>epsilon 
      Tcpu=NaN;
      NF=NaN;
      Itr=NaN;
      NG=NaN;
    else
      Tcpu=toc; 
      NG=norm(Fk);
    end
    finish=1;
%      Tcpu=toc;
end
% [Tcpu Itr NF norm(Fk) norm(dk)]
%Tcpu=cputime-t0;
