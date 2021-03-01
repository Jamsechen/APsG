%% **************************************************************************************************
                                      % Description
%% **************************************************************************************************
% Generate data for the model Y=XB+W, where  Y(n x q), X(n x m), B(m x q),
 %% input:
    % rho_X--correlation for X;
    % rho_E--correlation for W;
    % df--degree freedom of t-distribution;
    % et--error type; 
    % 1--normal distribution error; 2--t-distribution; error; 3--Laplce-distribution error;
    % cor--type of generating the covariance matrix of X
%% Output:
% B_true--the true coefficient matrix
% Y--XB_true + W

%% *************************************************************************************************
                                    % Main programming
%% *************************************************************************************************

function [X,SIGMA_X,Y_E,B_true]  = Generatedata(n,m,m0,q,rho_X,rho_E,df,et,cor)      

%% generate predictor
  if cor==1 
      for u=1:m
          for v=1:m
              SIGMA_X(u,v) = rho_X^(abs(u-v));
          end
      end
   end
   
   if cor==2
      SIGMA_X = rho_X*ones(m,m);
      for u=1:m
          SIGMA_X(u,u) = 1;
      end
   end
   
  MU_X = zeros(1,m);
  X    = mvnrnd(MU_X,SIGMA_X,n);
  % X = [ones(n,1) X];
      
%% generate true coefficient matrix
  B_true = [randn(m0,q); zeros(m-m0,q)];
  Y      = X*B_true; 
      
%% add noise
  MU_E    = zeros(1,q); 
  SIGMA_E = ones(q,q);
   for u=1:q
      for v=1:q
          SIGMA_E(u,v) = rho_E^(abs(u-v));
      end
   end
   if et==1
      ERROR = mvnrnd(MU_E,SIGMA_E,n);
            % multivariate Normal distribution
   end
   if et==2
      ERROR = 10*mvtrnd(SIGMA_E,df,n);
           % multivariate t-distribution
   end
   if et==3
      ERROR = generate_MVLaplce(q,n);
           % multivariate Laplce distribution
   end
  Y_E = Y + ERROR;
  