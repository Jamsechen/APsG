%% Output: C--B1��B2�ж�Ӧ�ж�Ϊ��ĸ���
%           IC--B1�е�������B2�ж�Ӧ�з���ĸ���
%           nozero--B1�з����еĸ���
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



