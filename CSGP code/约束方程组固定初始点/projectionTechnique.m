function get=projectionTechnique(n,x,nprob)
%功能：投影技术

switch nprob
    case 1 %exponential problem
        get = max(0,x);
    case 2 %modified logarithmic problem
        if sum(x) <= n && min(x) > -1
            get = x;
        else
            get=quadprog(speye(n),-2*x,ones(1,n),n,[],[],-ones(n,1),[]);
        end
    case 3 %nonsmooth problem
        if sum(x) <= n && min(x) >= 0
            get = x;
        else
            get=quadprog(speye(n),-2*x,ones(1,n),n,[],[],zeros(n,1),[]);
        end
    case 4 %strictly convex problem I
        get=max(0,x);
    case 5 %strictly convex problem II
        get=max(0,x);
    case 6 %nonsmooth problem II
        if sum(x) <= n && min(x) >= -1
            get = x;
        else
            get=quadprog(speye(n),-2*x,ones(1,n),n,[],[],-ones(n,1),[]);
        end
    case 7 %Discrete boundary value problem
        get=max(0,x);
    case 8 %Trigonometric problem
        get=max(-2,x);
    case 9 %Tridiagonal problem 
        get=max(0,x);
    case 10 %modified logarithmic problem
        get=max(-1,x);
    case 11 %nonsmooth problem
        get=max(0,x);
    case 12 %exponential problem
        get=max(0,x);
    case 13
        get=max(0,x);
    case 14 %nonsmooth problem II
         if sum(x) <= n && min(x) > -1
            get = x;
        else
            get=quadprog(speye(n),-2*x,ones(1,n),n,[],[],-ones(n,1),[]);
        end
    case 15 %modified tridiagonal problem
        get=max(0,x);
    case 16
        get=max(0,x);
    case 17 %strictly convex function II
        get=max(0,x);
    case 18
        get=max(0,x);
    case 19 %Trigexp function
        get=max(0,x);
    case 20 %modified problem
        get=max(0,x);
    case 21 %modified problem
        get=max(-3,x);
    case 22 %Penalty I
        get=max(-3,x);
    case 23  %nonsmooth problem
        if sum(x) <= n && min(x) >= 0
            get = x;
        else
            get=quadprog(speye(n),-2*x,ones(1,n),n,[],[],zeros(n,1),[]);
        end
    case 24 %modified trigonometric function
        get=max(-3,x);
    case 25 %tri exp function
        get=max(0,x);
    case 26  %logarithmicfunction
        get=max(0,x);
    case 27
         get=max(0,x);
    case 28 %nonsmooth problem 
         get=max(0,x);
    case 29
        if sum(x) <= n && min(x) >= 0
            get = x;
        else
            get=quadprog(speye(n),-2*x,ones(1,n),n,[],[],zeros(n,1),[]);
        end
    otherwise 
        disp('this problem is not exist!!! ')
end

end