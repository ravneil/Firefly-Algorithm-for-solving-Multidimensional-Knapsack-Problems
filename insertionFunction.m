function [S] = insertionFunction(x,l, v,n,limit)



Cost = l;
x_temp = x;    
y=x_temp;

% if v>1 && v<n
for i=1:v


i=randperm(n,2);%randperm randsample

 i1=i(1);
 i2=i(2);

if Cost>limit
  if y([i1])~=0 && y([i2])~=0
      y([i1])=0;
      y([i2])=0;
  end
    
    if y([i1])~=0  
  y([i1])=0;
  y([i2])=1;
    else
      y([i1])=1;
      y([i2])=0;
  end
 
%     y([i2])=0;
    
end
if Cost<limit
    
     if y([i1])~=1 && y([i2])~=1
      y([i1])=1;
      y([i2])=1;
     end
  
    if y([i1])~=1  
  y([i1])=1;
  y([i2])=0;
  else
      y([i1])=0;
      y([i2])=1;
  end
 
%     y([i2])=1;
end

    x_temp=y;
  

end



 S = x_temp;

