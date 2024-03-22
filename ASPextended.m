%% ASP extended model

clear all
n = 2; %number of jobs
n_u =7 ;
d = 7; %day length
c_w = 1;
c_i = 0;
c_l = 10; % cost of tardiness (overtime)
c_g = 0; % cost of earliness


% simulation
k = 25000; % number of scenarios
h = 2.*rand(n+n_u,k);
NS = binornd(ones(n+n_u,k),0.7);
% NS_h = h.*NS;
delay =  rand(n+n_u,k);
NS_delay_h = NS.*(h + delay);
p = ones(1,k).*1/k;
p_u= [0.7,0.5,0.4,0.3,0.2,0.1,0.05];

% cost vector
q = [c_w;c_i;c_l;c_g];

% initial x
x = zeros(n+n_u-1,1);

% results
[w_v2,x_v2,theta_v2,st_dev,CI] = dynamic_lshaped_single(h,q,n,n_u,k,p,p_u,x,d); 
