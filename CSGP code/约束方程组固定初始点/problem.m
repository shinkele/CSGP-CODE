function g = problem(n,x,nprob)
%功能：用于求解凸约束非线性方程组问题
%输入：n：自变量维数；x：自变量；npro：问题序号
%输出：g：函数值

    switch nprob
        case 1 %exponential problem
            g = zeros(n,1);
            g(1)=exp(x(1))-1;
            g(2:n)=exp(x(2:n))+x(2:n)-1;
        case 2 %strictly convex problem I
            g = exp(x)-1;
        case 3 %strictly convex problem II
            g = (1:n)'.*exp(x)./n-1;
        case 4 %nonsmooth problem 
            g = 2.*x-sin(abs(x));
        case 5 %nonsmooth problem II
            g = x - sin(abs(x-1));
        case 6 %modified logarithmic problem
            g = log(x+1)-x./n;
        case 7
            g = exp(x).^2+3.*sin(x).*cos(x)-1; 
        case 8 
            g = zeros(n,1);
            g(1)=2*x(1)+sin(x(1))-1;
            g(2:n-1)=2.*x(1:n-2)+2.*x(2:n-1)+sin(x(2:n-1))-1;
            g(n)=2*x(n)+sin(x(n))-1;
        case 9 %modified problem
            g = x - 2.*sin(abs(x-1));
        otherwise 
            disp('this problem is not exist!!! ')
    end
end

