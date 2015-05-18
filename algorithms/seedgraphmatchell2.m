function [ corr, corr_c ] = seedgraphmatchell2(A,B,s,convexStart)% ,alpha_type )
% Returns matching for the SGM problem using the SGM algorithm
% 
% [corr,corr_c] = seedgraphmatchell2(A,B,m,convexStar)
% 
% corr_c is the matching using only seed to nonseed data
% corr is the best matching using all data
%
% A,B are (s+n)x(s+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
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

if nargin < 4
    warning('convexStart not set, default is true');
    convexStart = true;
end


[totv,~]=size(A); % number of vertices
n=totv-s; % number of non-seeds
eyen=eye(n);
scale = 10000;
patience=25;
tol=.99;

%% Get seed->non-seed, non-seed->seed, and non-seed->non-seed adj matrices
A12=A(1:s,s+1:s+n);
A21=A(s+1:s+n,1:s);
A22=A(s+1:s+n,s+1:s+n);
B12=B(1:s,s+1:s+n);
B21=B(s+1:s+n,1:s);
B22=B(s+1:s+n,s+1:s+n);

%% Get the initial start
if( ~convexStart )
    % Use the baricenter
	P = ones(n)/n;
else
    % use the start from the convex relaxation
	[~,P]=relaxed_normAPPB_FW_seeds(A22,B22,s);
end

%% Matching just using seed to non-seed information
corr_c = lapjv(-(A21*B21'+A12'*B12), scale );%YiCaoHungarian( -(A21*B21'+A12'*B12) );%
corr_c=[ 1:s,  s+corr_c];


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
corr = lapjv(-P, scale);%YiCaoHungarian(-P);%


%init_vs_final_dissagreements = sum(ind0~=corr)

corr=[ 1:s,  s+corr];

