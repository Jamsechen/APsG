%% Output: C--B1与B2中对应行都为零的个数
%           IC--B1中的零行且B2中对应行非零的个数
%           nozero--B1中非零行的个数
function [C,IC] = finderror(B1,B2)
   m  = size(B1,1);
   C  = 0;
   IC = 0;
for i=1:m
    
    if (norm(B1(i,:),2)~=0)
       h1 = 1;
    else
       h1 = 0;
    end
   if (norm(B2(i,:))~=0)
       h2 = 1;
   else
       h2 = 0;
   end

if (h1==0)&(h2==0)
    C = C+1;
end
if (h1==0)&(h2~=0)
   IC = IC+1;  
end 
end



