%% ASP basic model 'brute-force'
clear all
n = 7; %number of jobs
d = 7; %day length
c_w = [ones(1,n-1).*1];
c_i = [ones(1,n-1).*0];
c_l = 10; % cost of tardiness (overtime)
c_g = 0; % cost of earliness


% simulation
k = 800; % number of scenarios
h = 2.*rand(n,k);
h(n,:) = h(n,:)-d;
p = ones(1,k).*1/k;

 
% % construction of T and W
W = zeros(n*k,2*n*k);

for j = 1:k
    W(1+(j-1)*n,1+(j-1)*n*2)=1;
for i = 2+(j-1)*n:(n-1)+(j-1)*n
    W(i,i-1+(j-1)*n)=-1;
    W(i,i+(j-1)*n)=1;
end
W(n*j,n-1+2*n*(j-1))=-1;
W(n*j,2*n-1+2*n*(j-1)) = 1;
W(n*j,2*n+2*n*(j-1)) = -1;
W(1+(j-1)*n:n*j-1,n+2*n*(j-1):(2*n-2)+2*n*(j-1))=-eye(n-1);
end
T = [eye(n-1);ones(1,(n-1)).*(-1)]; % T matrix
T = repmat(T,k,1);
W = [T,W];


% initialize x vector 
x = zeros(n-1,1);

% cost vector
 q = [c_w';c_i';c_l;c_g];
  
% vector of RHS constraints for each scenario
 b = reshape(h,k*n,1);

 % cost vector comprehensive of all scenarios
 q2 = [zeros(n-1,1);repmat(q,k,1).*p(1)];
 
 % all decision variables are positive
 lb = zeros(length(q2),1);
 
 % solution for the full lp problem
 sol = linprog(q2',[],[],W,b,lb,[]);
 x = sol(1:n-1);  % optimal x
 total_cost = q2'*sol; % total cost
 
 st_dev = std(q'*reshape(sol(n:n*2*k+(n-1),1),n*2,k)); % std of total costs
 CI = [total_cost-norminv(0.975)*st_dev/(sqrt(k)),total_cost+norminv(0.975)*st_dev/(sqrt(k))]; %confidence interval of total costs
 
 Expected_y = p*reshape(sol(n:n*2*k+(n-1),1),2*n,k)'; %vector of expected waiting times, idle times, lateness and earliness
 
 