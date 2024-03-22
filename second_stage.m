function second_stage = second_stage(n,k,x,d,h)
second_stage = zeros(n*2, k);
second_stage(1,:) = max(h(1,:)-x(1,1),0);
second_stage(n,:) = max(-h(1,:)+x(1,1),0);
for i = 2:n-1
    second_stage(i,:) = max(second_stage(i-1,:) + h(i,:) - x(i,1),0);
    second_stage(n-1+i,:) = max(-second_stage(i-1,:) - h(i,:) + x(i,1),0);
end

second_stage(2*n-1,:) = max(second_stage(n-1,:) + h(n,:) + sum(x) - d,0);
second_stage(2*n,:) = max(-second_stage(n-1,:) - h(n,:) - sum(x) + d,0);
end