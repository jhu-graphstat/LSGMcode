function [P,Pp]=relaxed_normAPPB_FW_seeds(A,B,seeds)
% this solves the problem min ||AD-DB||_F over
% doubly stochastic matrices D
% uses Frank Wolfe methodology, although faster convex quadratic program
% solvers are available, and could be used instead.
%
% A and B are the two graphs
% s is the number of seeds
% vertices 1,2,3,...,s are assumed to be the seeds vertices in each graph
% with the seeding taking vertex i to vertex i
verbose = 1;

AtA = A'*A;
BBt = B*B';

p=size(A,1);

f1 = @(P) norm(A*P-P*B,'fro')^2;

tol=5e-2;
tol2=1e-5;

P=ones(p)/(p-seeds);
P(1:seeds,1:seeds)=eye(seeds);

f=f1(P);
var=1;

%begin F-W iterates
while ~(abs(f)<tol) && (var > tol2)
    fold=f;

    grad = (AtA*P -A'*P*B - A*P*B' + P*BBt);
    
    grad(1:seeds,:)=0;
    grad(:,1:seeds)=0;
    
    corr=lapjv(grad(seeds+1:end,seeds+1:end),0.01);  
    Ps=perm2mat(corr);

    Ps =Ps';
    
    Ps1=eye(p);
    Ps1(seeds+1:end,seeds+1:end) = Ps;
    Ps=Ps1;
    
    C = A*(P-Ps) + (Ps-P)*B; 
    D = A*Ps-Ps*B;
    
    aq = trace(C*C');
    bq = trace(C*D'+D*C');
    aopt = -bq/(2*aq);

    Ps4 = aopt*P + (1-aopt)*Ps;
    
    f=f1(Ps4);
    P=Ps4;
    
    var=abs(f-fold);
    if (verbose) fprintf('f: %1.5f  var %.15f\n',f,var); end

end

    corr=lapjv(-P,0.01);
    Pp=perm2mat(corr);


end
