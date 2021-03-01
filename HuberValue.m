function hubervalue = HuberValue(X,huberC)
 n          = size(X,1);
 hubervalue = 0;
 for i=1:n
       C = norm(X(i,:),2);
    if C<=huberC
       hubervalue = hubervalue + C^2/2; 
    else
       hubervalue = hubervalue + huberC*(C-huberC/2);
    end
 end
