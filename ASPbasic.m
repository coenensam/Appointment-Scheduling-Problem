%% ASP basic model
clear all
tic
n = 9; %number of jobs
d = 7; %day length
c_w = 10;
c_i = 1;
c_l = 0; % cost of tardiness (overtime)
c_g = 0; % cost of earliness


% simulation
k = 2000; % number of scenarios
h = 2.*rand(n,k);
NS = binornd(ones(n,k),0.7);
delay =  2.*rand(n,k);
NS_h = h.*NS;
NS_delay_h = (h + delay);
p = ones(1,k).*1/k;



 
% % construction of T and W
W = zeros(n,2*n);

 W(1,1)=1;
for i = 2:(n-1)
    W(i,i-1)=-1;
    W(i,i)=1;
end
W(n,n-1)=-1;
W(n,2*n-1) = 1;
W(n,2*n) = -1;
W(1:n-1,n:(2*n-2))=-eye(n-1);

T = [eye(n-1);ones(1,(n-1)).*(-1)]; % T matrix



% initialize x vector 
x = zeros(n-1,1);

% cost vector
 q = [c_w;c_i;c_l;c_g];


[w_v,x_v,theta_v,st_dev,CI] = lshaped(T,NS_delay_h,q,n,k,p,x,d);  % single-cut

% 
% [w_v,x_v,theta_v,st_dev,CI] = lshaped_multi(T,h,q,n,k,p,x,d); % multi-cut


