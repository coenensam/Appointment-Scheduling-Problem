function T = T_construction (n)
T = [eye(n-1);ones(1,(n-1)).*(-1)]; % T matrix
end