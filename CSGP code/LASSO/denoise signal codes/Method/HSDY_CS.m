function [x,x_debias,objective,times,debias_start,mses,taus]= ...
    HSDY_CS(y,A,tau,varargin)
%
% Derivative-free HS-DY-type method for solving nonlinear equations and image restoration(2020)(HSDY_CS)
%
% This function solves the convex problem 
% arg min_x = 0.5*|| y - A x ||_2^2 + tau || x ||_1
% using the algorithm modified three-term conjugate gradient method, described in the following paper
%
% This code is to use the well-known code CG_DESCENT to solve \ell_1 norm 
% regularization least square problems.
% 
% -----------------------------------------------------------------------
% Copyright (2019): Jianghua Yin
% ----------------------------------------------------------------------
%
% 
% The first version of this code by Jianghua Yin, Nov. 29, 2019
% test for number of required parametres
if (nargin-length(varargin)) ~= 3
  error('Wrong number of required parameters');
end

% flag for initial x (can take any values except 0,1,2)
Initial_X_supplied = 3333;

% Set the defaults for the optional parameters
stopCriterion = 3;
tolA = 0.01;
tolD = 0.0001;
debias = 0;
maxiter = 10000;
maxiter_debias = 500;
miniter = 5;
miniter_debias = 5;
init = 0;
compute_mse = 0;
AT = 0;
verbose = 1;
continuation = 0;
cont_steps = -1;
firstTauFactorGiven = 0;

% Set the defaults for outputs that may not be computed
debias_start = 0;
x_debias = [];
mses = [];

% Read the optional parameters
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch upper(varargin{i})
     case 'STOPCRITERION'
       stopCriterion = varargin{i+1};
     case 'TOLERANCEA'       
       tolA = varargin{i+1};
     case 'TOLERANCED'
       tolD = varargin{i+1};
     case 'DEBIAS'
       debias = varargin{i+1};
     case 'MAXITERA'
       maxiter = varargin{i+1};
     case 'MAXITERD'
       maxiter_debias = varargin{i+1};
     case 'MINITERA'
       miniter = varargin{i+1};
     case 'MINITERD'
       miniter_debias = varargin{i+1};
     case 'INITIALIZATION'
       if prod(size(varargin{i+1})) > 1   % initial x supplied as array
	      init = Initial_X_supplied;      % flag to be used below
	      x = varargin{i+1};
       else 
	      init = varargin{i+1};
       end
     case 'MONOTONE'
       enforceMonotone = varargin{i+1};
     case 'CONTINUATION'
       continuation = varargin{i+1};  
     case 'CONTINUATIONSTEPS' 
       cont_steps = varargin{i+1};
     case 'FIRSTTAUFACTOR'
       firstTauFactor = varargin{i+1};
       firstTauFactorGiven = 1;
     case 'TRUE_X'
       compute_mse = 1;
       true = varargin{i+1};
     case 'ALPHAMIN'
       alphamin = varargin{i+1};
     case 'ALPHAMAX'
       alphamax = varargin{i+1};
     case 'AT'
       AT = varargin{i+1};
     case 'VERBOSE'
       verbose = varargin{i+1};
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%

if (sum(stopCriterion == [0 1 2 3 4 5])==0)
  error(['Unknown stopping criterion']);
end

% if A is a function handle, we have to check presence of AT,
if isa(A, 'function_handle') & ~isa(AT,'function_handle')
  error(['The function handle for transpose of A is missing']);
end 

% if A is a matrix, we find out dimensions of y and x,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~isa(A, 'function_handle')
  AT = @(x) (x'*A)'; %A'*x;
  A = @(x) A*x;
end
% from this point down, A and AT are always function handles.

% Precompute A'*y since it'll be used a lot
Aty = AT(y);

% Initialization
switch init
    case 0   % initialize at zero, using AT to find the size of x
       x = AT(zeros(size(y)));
    case 1   % initialize randomly, using AT to find the size of x
       x = randn(size(AT(zeros(size(y)))));
    case 2   % initialize x0 = A'*y
       x = Aty; 
    case Initial_X_supplied % initial x was given by user
       % initial x was given as a function argument; just check size
       if size(A(x)) ~= size(y)
          error(['Size of initial x is not compatible with A']); 
       end
    otherwise
       error(['Unknown ''Initialization'' option']);
end

% now check if tau is an array; if it is, it has to 
% have the same size as x
if prod(size(tau)) > 1
   try,
      dummy = x.*tau;
   catch,
      error(['Parameter tau has wrong dimensions; it should be scalar or size(x)']),
   end
end
     

% if the true x was given, check its size
if compute_mse & (size(true) ~= size(x))  
   error(['Initial x has incompatible size']); 
end

% if tau is scalar, we check its value; if it's large enough,
% the optimal solution is the zero vector
if prod(size(tau)) == 1
   aux = AT(y);
   max_tau = max(abs(aux(:)));
   if tau >= max_tau
      x = zeros(size(aux));
      if debias
         x_debias = x;
      end
      objective(1) = 0.5*(y(:)'*y(:));
      times(1) = 0;
      if compute_mse
          mses(1) = sum(true(:).^2);
      end
      return
   end 
end

% initialize u and v
u =  x.*(x >= 0);
v = -x.*(x <  0);

% define the indicator vector or matrix of nonzeros in x
nz_x = (x ~= 0.0);
num_nz_x = sum(nz_x(:));

% start the clock
t0 = tic;

% store given tau, because we're going to change it in the
% continuation procedure
final_tau = tau;

% store given stopping criterion and threshold, because we're going 
% to change them in the continuation procedure
final_stopCriterion = stopCriterion;
final_tolA = tolA;

% set continuation factors
if continuation&&(cont_steps > 1)
   % If tau is scalar, first check top see if the first factor is 
   % too large (i.e., large enough to make the first 
   % solution all zeros). If so, make it a little smaller than that.
   % Also set to that value as default
   if prod(size(tau)) == 1
      if (firstTauFactorGiven == 0)|(firstTauFactor*tau >= max_tau)
         firstTauFactor = 0.5*max_tau / tau;
         if verbose
             fprintf(1,'\n setting parameter FirstTauFactor\n')
         end
      end
   end
   cont_factors = 10.^[log10(firstTauFactor):...
                    log10(1/firstTauFactor)/(cont_steps-1):0];
end

if ~continuation
  cont_factors = 1;
  cont_steps = 1;
end
  
iter = 1;
if compute_mse
       mses(iter) = sum((x(:)-true(:)).^2);
end

keep_continuation = 1;
cont_loop = 1;
iter = 1;
taus = [];

% loop for continuation
while keep_continuation
    % Compute and store initial value of the objective function
    resid =  y - A(x);
    if cont_steps == -1
        gradq = AT(resid);
        tau = max(final_tau,0.2*max(abs(gradq)));
        if tau == final_tau
            stopCriterion = final_stopCriterion;
            tolA = final_tolA;
            keep_continuation = 0;                  % stop continuation
        else
            stopCriterion = 1;
            tolA = 1e-5;
        end
    else
        tau = final_tau * cont_factors(cont_loop);% 
        if cont_loop == cont_steps
            stopCriterion = final_stopCriterion;
            tolA = final_tolA;
            keep_continuation = 0;                 % 
        else
            stopCriterion = 1;
            tolA = 1e-5;
        end
    end
    taus = [taus tau];
    
    if verbose
        fprintf(1,'\nSetting tau = %0.5g\n',tau)
    end
    % if in first continuation iteration, compute and store
    % initial value of the objective function
    if cont_loop == 1
        alpha = 1.0;
        f = 0.5*(resid(:)'*resid(:)) + ...
             sum(tau(:).*u(:)) + sum(tau(:).*v(:));
        objective(1) = f;
        if compute_mse
            mses(1) = (x(:)-true(:))'*(x(:)-true(:));
        end
        if verbose
            fprintf(1,'Initial obj=%10.6e, alpha=%6.2e, nonzeros=%7d\n',...
                f,alpha,num_nz_x);
        end
    end
    % Compute the initial gradient and the useful 
    % quantity resid_base
    resid_base = y - resid;
    % control variable for the outer loop and iteration counter
    keep_going = 1;
    if verbose
      fprintf(1,'\nInitial obj=%10.6e, nonzeros=%7d\n',f,num_nz_x);
    end
    while keep_going 
      % compute gradient
      temp = AT(resid_base);
      term  =  temp - Aty;
      gradu =  term + tau; % Hz+c w.r.t. u
      gradv = -term + tau; % Hz+c w.r.t. v
      %
      Lu = min(u,gradu); % F(z) = min(z,Hz+c) w.r.t. u, z = [u v]'
      Lv = min(v,gradv); % F w.r.t. v
      NormF = Lu(:)'*Lu(:) + Lv(:)'*Lv(:); % compute ||F(z)||^2
      % NormF = sqrt(NormF); 
      % compute the search direction dk
      % which comes form the paper: 
      if (iter > 1)
%          % compute y_{k-1}
%          yku = Lu-old_Lu; 
%          ykv = Lv-old_Lv;
%          % compute y_{k-1}'*d_{k-1}
%          ykdk = yku(:)'*old_du(:)+ykv(:)'*old_dv(:);
%          % compute ||d_{k-1}|| and ||y_{k-1}||
%          NormD=sqrt(old_du(:)'*old_du(:)+old_dv(:)'*old_dv(:));
%          Normyk2=yku(:)'*yku(:)+ykv(:)'*ykv(:);
%          NormY=sqrt(Normyk2);
%          % compute w_{k}, F_{k}'*d_{k-1} and F_{k}'*y_{k-1}
%          wk=max(10*NormD*NormY,max(ykdk,Old_NormF)); % wk=max(0.2*NormD*NormY,max(rdk,Old_NormF));
%          Fkdk1 = Lu(:)'*old_du(:)+Lv(:)'*old_dv(:);
%          Fkyk1 = Lu(:)'*yku(:)+Lv(:)'*ykv(:);
%          % compute betak
%          Betak = Fkyk1/wk-Normyk2*Fkdk1/wk^2; 
%          % compute tk and lambdak
%          sk1u=u-old_u;
%          sk1v=v-old_v;
%          yksk1=yku(:)'*sk1u(:)+ykv(:)'*sk1v(:);
%          tk=min(0.3,max(0,1-yksk1/Normyk2)); % tk=min(0.5,max(0,1-yksk1/Normyk2));
%          lambdak = tk*Fkdk1/wk;
%          % computation of search direction vector
%          du = - Lu + Betak*old_du + lambdak*yku;
%          dv = - Lv + Betak*old_dv + lambdak*ykv;
         %=================================================================
         % compute sk, yk
         sku = u - old_u;
         skv = v - old_v;
         yku = Lu - old_Lu;
         ykv = Lv - old_Lv;
         % compute sk_bar, thetak
         skyk = [sku(:);skv(:)]'*[yku(:);ykv(:)];
         Normyk = norm([yku(:);ykv(:)]);
         max_term = max(0, -skyk/Normyk^2);
         sk_baru = sku + max_term * yku + yku;
         sk_barv = skv + max_term * ykv + ykv;
         
         yksk_bar = [sk_baru(:);sk_barv(:)]'*[yku(:);ykv(:)];
         thetak = Normyk^2/yksk_bar;
         % compute betak_mhs, betak_mdy
         dkyk = [old_du(:);old_dv(:)]'*[yku(:);ykv(:)];
         Normdk = norm([old_du(:);old_dv(:)]);
         tk = 1 + max(0,-dkyk/Normdk^2);
         wku = yku + tk * old_du;
         wkv = ykv + tk * old_dv;
         dkwk = [old_du(:);old_dv(:)]'*[wku(:);wkv(:)];
         
         Fk1yk = [Lu(:);Lv(:)]'*[yku(:);ykv(:)];
         betak_mhs = Fk1yk/dkwk;
         
         NormFk1 = norm([Lu(:);Lv(:)]);
         betak_mdy = NormFk1^2/dkwk;
         % compute betak_hsdy
         betak_hsdy = (1-thetak)*betak_mhs + thetak*betak_mdy;
                 
         % computation of search direction vector
         Fk1dk = [Lu(:);Lv(:)]'*[old_du(:);old_dv(:)];
         du = - (1+betak_hsdy*Fk1dk/NormFk1^2)*Lu + betak_hsdy*old_du;
         dv = - (1+betak_hsdy*Fk1dk/NormFk1^2)*Lv + betak_hsdy*old_dv;
         
      else
          du = - Lu;
          dv = - Lv;
      end
      % store 
      old_u = u; 
      old_v = v;
      old_du = du;
      old_dv = dv;
      Old_NormF = NormF;
%       old_gradu = gradu;
%       old_gradv = gradv;
      old_Lu = Lu;
      old_Lv = Lv;
      % compute dx
      dx = du-dv;
      % calculate useful matrix-vector product involving dx
      auv = A(dx); 
      Bdu = AT(auv);
      Bdv = -Bdu;
      % preparetion for line search 
      sigma = 0.0001;
      betas = 1;   %  betas = 10;
      Lu = min(u+betas*du, gradu+betas*Bdu); % F(z+betas*d) w.r.t. u where d=[du;dv];
      Lv = min(v+betas*dv, gradv+betas*Bdv); % F(z+betas*d) w.r.t. v;
%       NormFz = sqrt(Lu(:)'*Lu(:) + Lv(:)'*Lv(:));
      Luvduv = Lu(:)'*du(:) + Lv(:)'*dv(:); % F(z+betas*d)'*d
      dudv = du(:)'*du(:)+dv(:)'*dv(:); % ||d||^2
      % line search  process
      while - Luvduv < sigma*betas*dudv % - Luvduv < sigma*betas*max(0.001,min(0.8,NormFz))*dudv 
          betas = 0.8*betas;  % betas = 0.5*betas;
          Lu = min(u+betas*du,gradu+betas*Bdu);
          Lv = min(v+betas*dv,gradv+betas*Bdv);
%           NormFz = sqrt(Lu(:)'*Lu(:) + Lv(:)'*Lv(:));
          Luvduv = Lu(:)'*du(:) + Lv(:)'*dv(:);
      end
%       compute ||F(z+betas*d)||^2
%       NormFz2=Lu(:)'*Lu(:)+Lv(:)'*Lv(:);
      % compute the projection steplength
      lambda = -1.2*Luvduv*betas/(Lu(:)'*Lu(:)+Lv(:)'*Lv(:));  % lambda = - 1.6*Luvduv*betas/NormFz^2;
      %
      u = old_u - lambda * Lu;
      v = old_v - lambda * Lv;
%       if iter < 10
%       u = max(0,old_u - lambda * Lu);
%       v = max(0,old_v - lambda * Lv);
%       else
%       u = old_u - lambda * Lu;
%       v = old_v - lambda * Lv;
%       end
      uvmin = 0;% min(u,v);
      u = u - uvmin; 
      v = v - uvmin; 
      x = u - v;
      % calculate nonzero pattern and number of nonzeros (do this *always*)
      nz_x_prev = nz_x;
      nz_x = (x~=0.0);
      num_nz_x = sum(nz_x(:));
      % update residual and function
      ALuv = A(Lu-Lv);
      resid = y - resid_base + lambda*ALuv;
      prev_f = f;
      f = 0.5*(resid(:)'*resid(:)) +  sum(tau(:).*u(:)) + ...
          sum(tau(:).*v(:));
      % compute new alpha
      dd  = Lu(:)'*Lu(:) + Lv(:)'*Lv(:);  
      %
      resid_base = resid_base - lambda*ALuv; 
      % print out stuff
      if verbose
         fprintf(1,'It=%4d, obj=%9.5e, alpha=%6.2e, nz=%8d  ',...
             iter, f, alpha, num_nz_x);
      end
      % update iteration counts, store results and times
      iter = iter + 1;
      objective(iter) = f;
      times(iter) = toc(t0);

      if compute_mse
        err = true - x;
        mses(iter) = (err(:)'*err(:));
      end

      switch stopCriterion
          case 0,
              % compute the stopping criterion based on the change
              % of the number of non-zero components of the estimate
              num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
              if num_nz_x >= 1
                  criterionActiveSet = num_changes_active;
              else
                  criterionActiveSet = tolA / 2;
              end
              keep_going = (criterionActiveSet > tolA);
              if verbose
                  fprintf(1,'Delta n-zeros = %d (target = %e)\n',...
                      criterionActiveSet , tolA)
              end
          case 1,
              % compute the stopping criterion based on the relative
              % variation of the objective function.
              criterionObjective = abs(f-prev_f)/(prev_f);
              keep_going =  (criterionObjective > tolA);
              if verbose
                  fprintf(1,'Delta obj. = %e (target = %e)\n',...
                      criterionObjective , tolA)
              end
          case 2,
              % stopping criterion based on relative norm of step taken
              delta_x_criterion = norm(Lu(:)-Lv(:))/norm(x(:));
              keep_going = (delta_x_criterion > tolA);
              if verbose
                  fprintf(1,'Norm(delta x)/norm(x) = %e (target = %e)\n',...
                      delta_x_criterion,tolA)
              end
          case 3,
              % compute the "LCP" stopping criterion - again based on the previous
              % iterate. Make it "relative" to the norm of x.
              w = [ min(gradu(:), old_u(:)); min(gradv(:), old_v(:)) ];
              criterionLCP = norm(w(:), inf);
              criterionLCP = criterionLCP / ...
                  max([1.0e-6, norm(old_u(:),inf), norm(old_v(:),inf)]);
              keep_going = (criterionLCP > tolA);
              if verbose
                  fprintf(1,'LCP = %e (target = %e)\n',criterionLCP,tolA)
              end
          case 4,
              % continue if not yeat reached target value tolA
              keep_going = (f > tolA);
              if verbose
                  fprintf(1,'Objective = %e (target = %e)\n',f,tolA)
              end
          case 5,
            % stopping criterion based on relative norm of step taken
            delta_x_criterion = sqrt(dd)/sqrt(x(:)'*x(:));
            keep_going = (delta_x_criterion > tolA);
            if verbose
                fprintf(1,'Norm(delta x)/norm(x) = %e (target = %e)\n',...
                    delta_x_criterion,tolA)
            end
          otherwise,
              error(['Unknown stopping criterion']);
      end % end of the stopping criteria switch
      
      % take no less than miniter... 
      if iter<=miniter
	      keep_going = 1;
      elseif iter > maxiter %and no more than maxiter iterations  
	      keep_going = 0;
      end
      
   end % end of the main loop of keep_going

    % increment continuation loop counter
    cont_loop = cont_loop+1;
    
end % end of the continuation loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
   fprintf(1,'\nFinished the main algorithm!\nResults:\n')
   fprintf(1,'||A x - y ||_2^2 = %10.3e\n',resid(:)'*resid(:))
   fprintf(1,'||x||_1 = %10.3e\n',sum(abs(x(:))))
   fprintf(1,'Objective function = %10.3e\n',f);
   nz_x = (x~=0.0); num_nz_x = sum(nz_x(:));
   fprintf(1,'Number of non-zero components = %d\n',num_nz_x);
   fprintf(1,'CPU time so far = %10.3e\n', times(iter));
   fprintf(1,'\n');
end

% If the 'Debias' option is set to 1, we try to remove the bias from the l1
% penalty, by applying CG to the least-squares problem obtained by omitting
% the l1 term and fixing the zero coefficients at zero.

% do this only if the reduced linear least-squares problem is
% overdetermined, otherwise we are certainly applying CG to a problem with a
% singular Hessian

if (debias & (sum(x(:)~=0)~=0))
  
  if (num_nz_x > length(y(:)))
    if verbose
      fprintf(1,'\n')
      fprintf(1,'Debiasing requested, but not performed\n');
      fprintf(1,'There are too many nonzeros in x\n\n');
      fprintf(1,'nonzeros in x: %8d, length of y: %8d\n',...
	  num_nz_x, length(y(:)));
    end
  elseif (num_nz_x==0)
    if verbose
      fprintf(1,'\n')
      fprintf(1,'Debiasing requested, but not performed\n');
      fprintf(1,'x has no nonzeros\n\n');
    end
  else
    if verbose
      fprintf(1,'\n')
      fprintf(1,'Starting the debiasing phase...\n\n')
    end
    
    x_debias = x;
    zeroind = (x_debias~=0); 
    cont_debias_cg = 1;
    debias_start = iter;
    
    % calculate initial residual
    resid = A(x_debias);
    resid = resid-y;
    resid_prev = eps*ones(size(resid));
    
    rvec = AT(resid);
    
    % mask out the zeros
    rvec = rvec .* zeroind;
    rTr_cg = rvec(:)'*rvec(:);
    
    % set convergence threshold for the residual || RW x_debias - y ||_2
    tol_debias = tolD * (rvec(:)'*rvec(:));
    
    % initialize pvec
    pvec = -rvec;
    
    % main loop
    while cont_debias_cg
      
      % calculate A*p = Wt * Rt * R * W * pvec
      RWpvec = A(pvec);      
      Apvec = AT(RWpvec);
      
      % mask out the zero terms
      Apvec = Apvec .* zeroind;
      
      % calculate alpha for CG
      alpha_cg = rTr_cg / (pvec(:)'* Apvec(:));
      
      % take the step
      x_debias = x_debias + alpha_cg * pvec;
      resid = resid + alpha_cg * RWpvec;
      rvec  = rvec  + alpha_cg * Apvec;
      
      rTr_cg_plus = rvec(:)'*rvec(:);
      beta_cg = rTr_cg_plus / rTr_cg;
      pvec = -rvec + beta_cg * pvec;
      
      rTr_cg = rTr_cg_plus;
      
      iter = iter+1;
      
      objective(iter) = 0.5*(resid(:)'*resid(:)) + ...
	  sum(tau(:).*abs(x_debias(:)));
      times(iter) = toc(t0);
      
      if compute_mse
	err = true - x_debias;
	mses(iter) = (err(:)'*err(:));
      end
      
      % in the debiasing CG phase, always use convergence criterion
      % based on the residual (this is standard for CG)
      if verbose
	fprintf(1,' Iter = %5d, debias resid = %13.8e, convergence = %8.3e\n', ...
	    iter, resid(:)'*resid(:), rTr_cg / tol_debias);
      end
      cont_debias_cg = ...
	  (iter-debias_start <= miniter_debias )| ...
	  ((rTr_cg > tol_debias) & ...
	  (iter-debias_start <= maxiter_debias));
      
    end
    if verbose
      fprintf(1,'\nFinished the debiasing phase!\nResults:\n')
      fprintf(1,'||A x - y ||_2^2 = %10.3e\n',resid(:)'*resid(:))
      fprintf(1,'||x||_1 = %10.3e\n',sum(abs(x(:))))
      fprintf(1,'Objective function = %10.3e\n',f);
      nz = (x_debias~=0.0);
      fprintf(1,'Number of non-zero components = %d\n',sum(nz(:)));
      fprintf(1,'CPU time so far = %10.3e\n', times(iter));
      fprintf(1,'\n');
    end
  end
  
  if compute_mse
    mses = mses/length(true(:));
  end
  
end

  


