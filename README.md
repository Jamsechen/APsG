# This is a matlab code for solving row Elastic-Net regularized multivariate Huber regression model 
 
 Model: 
         min_B {sum(h_t(Y_{i.}-X_{i.}B) + mu1* sum(||B_{i.}||) +(mu2/2)* ||B||_F^2}
 
 where h() is the vector form Huber function and B_{i.} is the ith row of B, and h_t(x)=||x||^2/2, when ||x||<=t; h_t(x)=t(||x||-t/2), when ||x||>=t.
 
 The main code is "APsG_Huber.m", which is in the form 
      
      [B_new,Supp_B_new,Time] = APsG_Huber(X,Y,B0,mu1,mu2,huberC,para)

Input:
      
      X--Prediction matrix;
      
      Y--Response matrix;
      
      B0--Start point;
      
      mu1--tuning parameter for l1/l2 term;
      
      mu2--tuning parameter for Frobenius term;
      
      huberC--Huber constant;
      
      para--other parameters 
          
            para.maxiter--the maximum number of iterations;
          
            para.eta--linesearch parameter
          
            para.tol--tolerance
          
            para.linesearch--control linesearch (1--linesearch; 0--no linesearch)

Output:
      
      B_new--the estimator
      
      Supp_B_new--The index of nonzero row of B_new
      
      Time--CPU time

Example:

In the file "Runexample.m", we present an example for run our code.
