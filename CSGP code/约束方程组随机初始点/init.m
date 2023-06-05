function x0 = init(n, type)
% 功能：产生初始点
% 输入：自变量维数 n；初始点类型 type
% 输出：初始点 x0

    switch type
        case 1
            x0 = ones(n,1); % (1, 1, ..., 1)
        case 2
            x0 = (1./logspace(log10(3),log10(3.^n),n))'; %(1/3, 1/3^2, 1/3^3, ..., 1/3^n)
        case 3
            x0 = (1./logspace(log10(2),log10(2.^n),n))'; %(1/2, 1/2^2, 1/2^3, ..., 1/2^n)
        case 4
            x0 = ((linspace(1,n,n)-1)./n)'; % (0, 1/n, 2/n, ..., (n-1)/n);
        case 5
            x0 = (1./linspace(1,n,n))'; %(1, 1/2, 1/3, ..., 1/n)
        case 6 
            x0 = (linspace(1,n,n)./n)'; %(1/n, 2/n, ..., n/n)
        case 7
            x0 = 1-(linspace(1,n,n)./n)'; %(1-1/n, 1-2/n, ..., 1-n/n)
        otherwise
            disp('没有该初始点！！！');
    end

end