%% min_B { sum[h(Y_{i.}-X_{i.}B] + mu1*sum[||B_{i.}||_2] + (mu2/2)*||B||_F^2 }
%%      where h() is the vector form Huber function and B_{i.} is the ith row of B.
%% h_t(x)=||x||^2/2, when ||x||<=t; h_t(x)=t(||x||-t/2), when ||x||>=t.

function [B_new,Supp_B_new,Time] = APsG_Huber(X,Y,B0,mu1,mu2,huberC,para)
  n          = size(Y,1);
  m          = size(B0,1);
  Maxiter    = para.maxiter;
  eta        = para.eta;
  tol        = para.tol;
  linesearch = para.linesearch; 
  L_f        = eigs(X'*X,1,'LM');
  L_min      = 1e-2*L_f;

  Z_old = inv(mu2/n*eye(m)+X'*X)*X'*Y;
  t_old = 1;
  L_old = L_f;
  B_old = B0;

  time_start= clock;
for i=1:Maxiter

%% no linesearch
if linesearch==0
   L         = L_old;
   GradZ_old = Gradhuber(X,Y,Z_old,huberC); 
   G         = (L*Z_old - GradZ_old)/(L+mu2); 
   [B_new,Supp_B_new] = Shrinkage_Block(G,mu1/(L+mu2));
   B_new     = (1+mu2/n)*B_new; %If the predictors are standardized, change (1+mu2/n) to (1+mu2)
   L_new     = L;
end

%% linesearch

if linesearch==1
   L0 = max(eta*L_old,L_min);
for Lmax=1:50
      L            = L0;
      GradZ_old    = Gradhuber(X,Y,Z_old,huberC); 
      G            = (L*Z_old - GradZ_old)/(L+mu2); 
[B_temp,Supp_B_temp] = Shrinkage_Block(G,mu1/(L+mu2));
      B_temp       = (1+mu2/n)*B_temp; %If the predictors are standardized, change (1 + mu2/n) to (1+mu2)
      normB_temp2  = 0;
      normB_tempF2 = 0;
for j=1:m
    normj        = norm(B_temp(j,:),2);
    normB_temp2  = normB_temp2 + normj;
    normB_tempF2 = normB_tempF2 + normj^2;
end
    F_temp = HuberValue(Y-X*B_temp,huberC) + mu1*normB_temp2 + mu2*normB_tempF2/2;
    Q_temp = HuberValue(Y-X*Z_old,huberC) + trace(GradZ_old'*(B_temp-Z_old))... 
                           + L*norm(B_temp-Z_old,'fro')^2/2 + mu1*normB_temp2 + mu2*normB_tempF2/2;
if (F_temp <= Q_temp)|(Lmax==50)
    L_new      = min(L,L_f);
    B_new      = B_temp; 
    Supp_B_new = Supp_B_temp; 
    break; 
else
    L0 = min(L/eta,L_f);
end
end
end
%%
   t_new = (1+sqrt(1+4*t_old^2))/2;
   c     = (t_old-1)/(t_new);
   Z_new = (1+c)*B_new - c*B_old;

%% Check stopping criterion 
   GradB_new = Gradhuber(X,Y,B_new,huberC);
   GradZ_new = Gradhuber(X,Y,Z_new,huberC);
   R = L_new*(Z_new-B_new)-(GradZ_new - GradB_new);
   Cr1 = norm(R,'fro')/(L_new*max(1,norm(B_new,'fro')));
   Cr2 = abs(norm(Y-X*B_new,'fro')-norm(Y-X*B_old,'fro'))/max(1,norm(Y,'fro'));
if (Cr1<= tol)&(Cr2<5*tol)
   break;
else
    B_old = B_new;
    Z_old = Z_new;
    t_old = t_new;
    L_old = L_new;
end
end
   Time = etime(clock,time_start);
