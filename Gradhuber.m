%% The gradient of Huber(Y-XB) with respect to B

function Grad = Gradhuber(X,Y,B,alpha) 
    n          = size(Y,1);
    Res        = Y-X*B;
    Grad_Temp  = Res;
for i=1:n
    C = alpha/norm(Res(i,:),2);
    if C<1
       Grad_Temp(i,:) = C* Res(i,:);
    end
end
    Grad = -X'*Grad_Temp;