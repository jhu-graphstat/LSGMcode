% Convert a symmetrical adjacency matrix to a Pajek network file
%
% Input
%   - adj : adjacency matrix
%   - name: file name (without extension)
%   - dir : destination directory
%
% Author: Erwan Le Martelot
% Date: 31/05/11

function [] = adj2pajek(adj, name, dir)

    dst_file = [name,'.net'];
    if nargin == 3
        dst_file = [dir, '/', dst_file];
    end

    fid = fopen(dst_file, 'w');
    
    fprintf(fid, '*Vertices %d\n', length(adj));
    
    for i=1:length(adj)
        fprintf(fid, ' %d "v%d"\n', i, i);
    end
    
    %fprintf(fid, '*Arcs\n');
    fprintf(fid, '*Edges\n');
    if ~issparse(adj)
        for i=1:length(adj)
            for j=i+1:length(adj)
                if adj(i,j) > 0
                    fprintf(fid, ' %d %d %f\n', i, j, adj(i,j));
                end
            end
        end
    else
        [i,j] = find(adj);
        for ind = 1:length(i)
            fprintf(fid, ' %d %d %f\n', i(ind), j(ind), full(adj(i(ind),j(ind))));
        end
    end
    fclose(fid);

end
