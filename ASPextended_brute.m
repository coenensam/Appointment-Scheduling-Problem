%% ASP extended model 'brute-force'
clear all
tic
n = 3; %number of jobs
d = 7; %day length
c_w = 1;
c_i = 0;
c_l = 10; % cost of tardiness (overtime)
c_g = 0; % cost of earliness


% simulation
n_u = 5;
k = 300; % number of scenarios
h = 2.*rand(n+n_u,k);
% h(n,:) = h(n,:)-d;
p = ones(1,k).*1/k;
p_u= [0.7,0.4,0.3,0.1,0.1,0.1,0.1,0.1,0.1];


p_stage = zeros(1,n_u+1);
p_stage(1) = 1 - p_u(1);
for i = 2:n_u
    p_stage(i) = prod(p_u(1:i-1))*(1-p_u(i));
end
p_stage(end) = prod(p_u);

 
% % construction of T and W
W1 = zeros([n:1:n+n_u]*ones(n_u+1,1).*k,2*[n:1:n+n_u]*ones(n_u+1,1).*k);
a = 0;
b = 0;
for r = n:n+n_u
    W = zeros(r*k,2*r*k);
for j = 1:k
    W(1+(j-1)*r,1+(j-1)*r*2)=1;
for i = 2+(j-1)*r:(r-1)+(j-1)*r
    W(i,i-1+(j-1)*r)=-1;
    W(i,i+(j-1)*r)=1;
end
W(r*j,r-1+2*r*(j-1))=-1;
W(r*j,2*r-1+2*r*(j-1)) = 1;
W(r*j,2*r+2*r*(j-1)) = -1;
W(1+(j-1)*r:r*j-1,r+2*r*(j-1):(2*r-2)+2*r*(j-1))=-eye(r-1);
end
W1(1+a:a+r*k,1+a*2:a*2+r*k*2) = W;
a = a+r*k;
end
T = [];
for i = n:n+n_u
Ttemp = [eye(i-1);ones(1,(i-1)).*(-1)]; % T matrix
Ttemp = repmat(Ttemp,k,1);
Ttemp = [Ttemp,zeros(i*k,n+n_u-i)];
T = [T;Ttemp];
end
W = [T,W1];


% initialize x vector 
x = zeros(n-1,1);
q = zeros(1,([n-1:1:n+n_u-1]*ones(n_u+1,1).*2+2*(n_u+1)).*k);
% cost vector
a = 0;
b = 0;
 for i = n:n+n_u
     a = b+((i-1)*2+2)*k;
     q(1,1+b:a) = repmat(p_stage(i-n+1).*p(1).*[c_w*ones(1,i-1),c_i*ones(1,i-1),c_l,c_g],1,k);
     b = a;
 end
 
  
% vector of RHS constraints for each scenario
b = [];
 for i = n:n+n_u
     btemp = h(1:i,:);
     btemp(end,:) = btemp(end,:) - d;
     btemp = reshape(btemp,i*k,1);
     b = [b;btemp];
 end

 % cost vector comprehensive of all scenarios
q2 = [zeros(1,n+n_u-1),q];
 
 % all decision variables are positive
 lb = zeros(length(q2),1);
 
 % solution for the full lp problem
 sol = linprog(q2,[],[],W,b,lb,[]);
 x = sol(1:n+n_u-1);  % optimal x
 total_cost = q2*sol; % total cost

 
 a=0;
 b=0;
 tot_cost = zeros(k,1);
 for i = n:n+n_u
     a = b+((i-1)*2+2)*k;
tot_cost = tot_cost + sum(reshape(q(1,1+b:a).*sol(n+n_u+b:n+n_u+a-1,1)'./p(1),k,((i-1)*2+2)),2);
b = a;
 end
 st_dev = std(tot_cost);
 CI = [total_cost-norminv(0.975)*st_dev/(sqrt(k)),total_cost+norminv(0.975)*st_dev/(sqrt(k))]; %confidence interval of total costs
 
