function [P,Pp]=relaxed_normAPPB_FW_seeds(A,B,seeds,roundGrad,varargin)

verbose = 0;

AtA = A'*A;
BBt = B*B';

p=size(A,1);

f1 = @(P) norm(A*P-P*B,'fro')^2;

tol=5e-2;
tol2=1e-5;
scale = 100;

% adds the ability to have seeds=0 in the function
if ( (numel(seeds)==1) & (seeds ~=0) )
    warning('Defaulting seeds to be the first %i',seeds);
    seeds = 1:seeds;
end
if ( (numel(seeds)==1) & (seeds ==0) )
    seeds=[];
end

nonSeeds = ~ismember(1:p,seeds);

nSeed = numel(seeds);
nnSeed = p-nSeed;

P = eye(p);
P(nonSeeds,nonSeeds)=ones(p-nSeed)/(p-nSeed);

if nargin < 4
    roundGrad = true;
end

if ~isempty(varargin)
    if (size(varargin{1},1) > 1)
        P=varargin{1}; 
    else
        tol2=varargin{1};    
    end
end

f=f1(P);
var=1;

iter = 0;
Pdiff = Inf;

while (f > tol) && (var > tol2) && (Pdiff>tol) 
    iter = iter+1;
    fold=f;

    grad = AtA*P -A'*P*B - A*P*B' + P*BBt;
    grad = grad + min(min(grad));
    grad(seeds,:)=0;
    grad(:,seeds)=0;
    if roundGrad
        corr=lapjv(round(nnSeed*grad(nonSeeds,nonSeeds)));
    else
        corr=lapjv(grad(nonSeeds,nonSeeds));
    end
    Ps=perm2mat(corr);

    Ps =Ps';
    
    Ps1=eye(p);
    Ps1(nonSeeds,nonSeeds) = Ps;
    Ps=Ps1;
    
    C = A*(P-Ps) + (Ps-P)*B; 
    D = A*Ps-Ps*B;
    
    aq = trace(C*C');
    bq = trace(C*D'+D*C');
    aopt = -bq/(2*aq);

    Ps4 = aopt*P + (1-aopt)*Ps;
    
    f=f1(Ps4);
    Pdiff = sum(sum(abs(P-Ps4)));
    P=Ps4;
    
    
    
    var=abs(f-fold);
    if (verbose) 
        fprintf('f: %1.5f  var %.15f\n',f,var);
    end

end

    corr=lapjv(-P,0.01);
    Pp=perm2mat(corr);


end
