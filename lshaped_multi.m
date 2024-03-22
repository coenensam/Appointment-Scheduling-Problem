%% multi-cut L-shaped method
function [w_v,x_v,theta_v,st_dev,CI] = lshaped_multi(T,h2,q,n,k,p,x,d) 
h = h2;
h(end,:) = h(end,:) - d;
c_i = q(2,1);
c_w = q(1,1);
c_l = q(3,1);
Dual_multip = dual_recursion(h2,c_i,c_w,c_l,n,k,x,d);

b=-sum(Dual_multip.*h)';
A_temp = Dual_multip'*T;
A = [-A_temp,-eye(k)];
lbx = zeros(n-1+k,1);
condition = 0;

while (condition == 0)
Master_c = [zeros(1,n-1),p];
Master_sol = linprog(Master_c,A,b,[],[],lbx,[]);

theta_v = p*Master_sol(n:k+n-1,1); % theta from the master is the lower bound
% objective primal >= objective dual, this relation yields the optimality cuts 
x = Master_sol(1:n-1,1);


Dual_multip = dual_recursion(h2,c_i,c_w,c_l,n,k,x,d);

Dual_c=h - T*x;
w_v = p*(sum(Dual_c.*Dual_multip))';

 
% check whether upper bound is greater than lower bound
  if((w_v-theta_v)>0)   
        A_temp = Dual_multip'*T; % if so, add optimality cuts
        A = [A;-A_temp,-eye(k)];
        b = [b;-(sum(Dual_multip.*h))'];
       
 else
     condition = 1;
          st_dev = std((sum(Dual_c.*Dual_multip))); % std of total costs
 CI = [w_v-norminv(0.975)*st_dev/(sqrt(k)),w_v+norminv(0.975)*st_dev/(sqrt(k))]; %confidence interval of total costs
 
 end
 
end
x_v = x;
end