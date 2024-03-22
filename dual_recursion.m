function f = dual_recursion(h,c_i,c_w,c_l,n,k,x,d) 
second_stage = zeros(n*2, k);
second_stage(1,:) = max(h(1,:)-x(1,1),0);
second_stage(n,:) = max(-h(1,:)+x(1,1),0);
for i = 2:n-1
    second_stage(i,:) = max(second_stage(i-1,:) + h(i,:) - x(i,1),0);
    second_stage(n-1+i,:) = max(-second_stage(i-1,:) - h(i,:) + x(i,1),0);
end

second_stage(2*n-1,:) = max(second_stage(n-1,:) + h(n,:) + sum(x) - d,0);
second_stage(2*n,:) = max(-second_stage(n-1,:) - h(n,:) - sum(x) + d,0);

dual_multip = zeros(n,k);
temp = second_stage(2*n-1,:);
temp(temp(1,:)>0)=c_l;
dual_multip(n,:) = temp;

for i = -(n-1):-1
    temp1 = c_w + dual_multip(-i+1,:);
    temp2 = (second_stage(-i,:)==0);
    temp1(temp2) = -c_i;
    dual_multip(-i,:) = temp1;
    
end
f = dual_multip;
end


    
    