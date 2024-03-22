 %% L-shaped single-cut for dynamic scheduling
function [w_v,x,theta_v,st_dev,CI] = dynamic_lshaped_single(h2,q,n,n_u,k,p,p_u,x,d) % single cut
cuts = cell(2,n_u+1);
w_v = zeros(1,n_u+1);
A = [];
b = [];
c_w = q(1);
c_i = q(2);
c_l = q(3);


p_stage = zeros(1,n_u+1);
p_stage(1) = 1 - p_u(1);
for i = 2:n_u
    p_stage(i) = prod(p_u(1:i-1))*(1-p_u(i));
end
p_stage(end) = prod(p_u);



condition = 0;
while(condition == 0)
Atemp = 0;
btemp = 0;
E = cell(1,n_u+1);
h3 = [];
dual_sub = [];
TT = [];
for i = n:n_u+n
T = [T_construction(i),zeros(i,n_u+n-i)];
dual_sub = [dual_sub;p_stage(i-n+1).*dual_recursion(h2(1:i,:),c_i,c_w,c_l,i,k,x(1:i-1,1),d)];
h3 = [h3;h2(1:i,:)];
h3(end,:) = h3(end,:) - d;
TT = [TT;T];
end 

btemp = p*(sum(dual_sub.*h3))';
Atemp = p*dual_sub'*TT;

Dual_c=h3 - TT*x;
w_v= p*(sum(Dual_c.*dual_sub))';

Atemp = [-Atemp,-1];
A = [A;Atemp];
b = [b;-btemp];

lbx = zeros(n+n_u,1);
Master_c = [zeros(1,n+n_u-1),1];  
Master_sol = linprog(Master_c,A,b,[],[],lbx,[]);
 
theta_v = Master_sol(n+n_u:end,1);
x = Master_sol(1:n+n_u-1,1);


if((w_v - theta_v)/theta_v<0.00001)
    condition = 1;
      st_dev = std((sum(Dual_c.*dual_sub))); % std of total costs
 CI = [w_v-norminv(0.975)*st_dev/(sqrt(k)),w_v+norminv(0.975)*st_dev/(sqrt(k))]; %confidence interval of total costs
 
end
end
end