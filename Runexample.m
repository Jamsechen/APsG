  clear all; 
  n  = 100;     m = ceil(3*n^(1/3)); 
  m0 = 7;       q = 5;
  
  rho_E = 0.5;   cor = 1;     
  df    = 2; %the degree freedom of t-distribution error               
        
  et    = 1; % 1--normal error; 2--t-distribution error;           
  rho_X = 0.9; 
  
  para.linesearch = 1;    % 1--linesearch; 0--no linesearch;
  para.maxiter    = 100;  %maximum number of iterations 
  para.eta        = 0.8;                   
  para.tol        = 1e-4;    
                    
  huberC = 10:-1:1;
  lhuber = length(huberC);
  
  mu_max = 200;   delta = 5*10^(-3);   mu_N = 100; 
  mu     = mu_max*delta.^([0:1:(mu_N-1)]/(mu_N-1));
  lmu    = length(mu);
  alpha  = 0.001;
  mu1    = (1-alpha)*mu;  
  mu2    = alpha*mu;
  
    N = 100;
for l=1:N
    [X,Sigma_X,Y_E,B_true] = Generatedata(n,m,m0,q,rho_X,rho_E,df,et,cor);
    B0           = zeros(size(B_true)); 
    normB_true2  = 0;
    normB_trueF2 = 0;
for j=1:m
    norm2        = norm(B_true(j,:),2);
    normB_true2  = normB_true2 + norm2;
    normB_trueF2 = normB_trueF2 + norm2^2;
end
       

    B_old = B0;
for k=1:lmu
    for i=1:lhuber
        [B,Supp_B,Cpu] = APsG_Huber(X,Y_E,B_old,mu1(k),mu2(k),huberC(i),para);
        Supp(k,i)      = length(Supp_B);
        Time(k,i)      = Cpu;
        [C_temp(k,i),IC_temp(k,i)] = finderror(B,B_true);
        MSE_temp(k,i) = trace((B-B_true)'*Sigma_X*(B-B_true));
        AIC(k,i)      = log(HuberValue(Y_E-X*B,huberC(i))/n) + Supp(k,i)*log(n)/n;
        B_old         = B;
    end
end

%% Select tunning parameter by Adjusted AIC
   [Rnz,Cnz] = find(Supp==0);
for k=1:length(Rnz)
    AIC(Rnz(k),Cnz(k)) = 1e+3;
end
    min_AIC     = min(min(AIC));
    [Row,Col]   = find(AIC==min_AIC);
    MSE_APSG(l) = MSE_temp(Row(1),Col(1));
    C_APSG(l)   = C_temp(Row(1),Col(1));
    IC_APSG(l)  = IC_temp(Row(1),Col(1));
    CPU_APSG(l) = Time(Row(1),Col(1));
end
  
%% Results after 100 repeat
   APSG_MSE = mean(MSE_APSG);      APSG_MSE_std = std(MSE_APSG);
   APSG_C   = mean(C_APSG);        APSG_C_std   = std(C_APSG);
   APSG_IC  = mean(IC_APSG);       APSG_IC_std  = std(IC_APSG);
   APSG_CPU = mean(CPU_APSG);      APSG_CPU_std = std(CPU_APSG);
%% Display the final result
 Result = {'Measurement', 'Estimation MSE', 'std(Estimation MSE)',      'C',    'std(C)',      'IC',    'std(IC)',     'CPU',    'std(CPU)';...
                    '--',         APSG_MSE,          APSG_MSE_std,   APSG_C,  APSG_C_std,   APSG_IC,  APSG_IC_std,  APSG_CPU,  APSG_CPU_std;};
 Result
