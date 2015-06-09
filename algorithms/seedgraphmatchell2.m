function [ corr, corr_c ] = seedgraphmatchell2(A,B,s,topK, start)% ,alpha_type )
% Returns matching for the SGM problem using the SGM algorithm
% 
% [corr,corr_c] = seedgraphmatchell2(A,B,m,start)
% 
% corr_c is the matching using only seed to nonseed data
% corr is the best matching using all data
%
% A,B are (s+n)x(s+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
% s is the number of seeds
% start is the initialization parameter
% It is assumed that the first s vertices of A's graph
% correspond respectively to the first m vertices of B's graph,
% corr gives the vertex correspondences  
% For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30 
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(400)]; B=round(rand(405,405));A=B(v,v);
% >> [corr,P] = seedgraphmatchell2( A,B,5 ) ; [v; corr]
% Extends Donniell's code
% (Extends Vogelstein, Conroy et al method for nonseed graphmatch to seed)

if nargin < 5
    warning('Input variable start not set, default is "convex"');
    start = 'convex';
end
if numel(seeds) == 1
    warning('Defaulting to seeds being the first %i vertices',seeds)
    seeds = 1:seeds;
end

s = numel(seeds);


[totv,~]=size(A); % number of vertices
n=totv-s; % number of non-seeds
eyen=eye(n); 
scale = 10000;
patience=25;
tol=.99;

nonSeeds = true(1,totv);
nonSeeds(seeds) = false;

%% Get seed->non-seed, non-seed->seed, and non-seed->non-seed adj matrices
A12=A(seeds,nonSeeds);
A21=A(nonSeeds,seeds);
A22=A(nonSeeds,nonSeeds);
B12=B(seeds,nonSeeds);
B21=B(nonSeeds,seeds);
B22=B(nonSeeds,nonSeeds);

%% Get the initial start
if( strcmp(start,'bari' ))
    % Use the baricenter
	P = ones(n)/n;
elseif( strcmp(start,'convex'))
    fprintf('Start Conv Relax\n')
    % use the start from the convex relaxation
	[~,P]=relaxed_normAPPB_FW_seeds(A22,B22,s);
    fprintf('Done Conv Relax\n')
else
    P = start;
end

%% Matching just using seed to non-seed information
if s > 0
    Ptotv = zeros(totv);
    Ptotv(nonSeeds,nonSeeds) = A21*B21'+A12'*B12;
    Ptotv(seeds,seeds) = eye(numel(seeds));
    corr_c = lapjv(-Ptotv, scale );%YiCaoHungarian( -Ptotv );%

else
    corr_c = NaN;
end

%% The main algorithm
toggle=1;
iter=0;
while (toggle==1)&&(iter<patience)
    iter=iter+1;
    % Compute the gradient of the objective function
    Grad=A22*P*B22'+A22'*P*B22+A21*B21'+A12'*B12;
    % Find the LAP solution for the negative gradient
    ind = lapjv( -Grad, scale );%YiCaoHungarian(-Grad);%
    T=eyen(ind,:);
    
%    if (alpha_type == 2)
%	    alpha = 2/(2+iter);
%	    P = (1-alpha)*P+(alpha)*T;
%    else

        % Compute some temporary quantities
		c=trace(A22'*P*B22*P');
		d=trace(A22'*T*B22*P')+trace(A22'*P*B22*T');
		e=trace(A22'*T*B22*T');
		u=trace(P'*A21*B21'+P'*A12'*B12);
		v=trace(T'*A21*B21'+T'*A12'*B12);
        
        % Determine step size
		alpha=-(d-2*e+u-v)/(2*(c-d+e));
		f0=0;
		f1=c-e+u-v;
		falpha=(c-d+e)*alpha^2+(d-2*e+u-v)*alpha;
        
        % Take the right step
		if (alpha<tol)&&(alpha>0)&&(falpha>f0)&&(falpha>f1)
		    P=alpha*P+(1-alpha)*T;
		elseif (f0>f1)
		    P=T;
        else
            % Stop now
		    toggle=0;
		end
%	end
end

%% Get the final correspondence
Ptotv = zeros(totv);
Ptotv(nonSeeds,nonSeeds) = P;
Ptotv(seeds,seeds) = eye(numel(seeds));
corr = lapjv(-Ptotv, scale);%YiCaoHungarian(-P);%


% %init_vs_final_dissagreements = sum(ind0~=corr)
% corr = 1:totV;
% corr(seeds) = seeds;
% corr(nonSeeds) = corrNS;

