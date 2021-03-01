# This is a matlab code for solving multivariate Huber regression model 
 min_B { sum[h(Y_{i.}-X_{i.}B] + mu1*sum[||B_{i.}||_2] +(mu2/2)*||B||_F^2}
 where h() is the vector form Huber function and B_{i.} is the ith row of B.
 h_t(x)=||x||^2/2, when ||x||<=t; h_t(x)=t(||x||-t/2), when ||x||>=t.