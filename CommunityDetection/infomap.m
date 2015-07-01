% InfoMap method from (Roswall and Bergstrom, 2008)
%
% Input
%   - adj: adjacency matrix
%
% Output
%   - com: communities (listed for each node)
%
% Author: Erwan Le Martelot
% Date: 16/06/11

function com = infomap(adj,trials, lambda)

    % Set the path and command line name 
    dir_path = '~/Documents/LSGM_code/CommunityDetection/clustering/scripts/';
    command = 'infomap';

    command = [dir_path,command];
    command = [command,' ',num2str(randi(99999999)),' adj.net',' ',num2str(trials),' ',num2str(lambda)];
    %command = [command, ' -d',' adj.net',' ./'];
    
    % Get edges list and output to file
    addpath('CommunityDetection')
    adj2pajek(adj,'adj','.');
 
    % Call community detection algorithm
    tic;
    system(command);
    toc;
    
    % Load data and create data structure
    % for the latest infomap 
%     fid = fopen('adj.tree','rt');
%     results = textscan(fid,'%s %f %q %f','Delimiter',' ','CommentStyle','#');
%     fclose(fid);
%     class = results{1};
%     Q = zeros(length(class),1);
%     for i = 1:length(class)
%         Q(i) = str2num(class{i}(1)); % first level class
%     end
%     nodeID = results{4}; 
%     [tmp,I] = sort(nodeID);
%     com = Q(I);
    
    % for the older infomap
    fid = fopen('adj.clu','rt');
    com = textscan(fid,'%d','CommentStyle','*');    
    fclose(fid);
    com = com{1};

    % Delete files
    %delete('adj.net');
    delete('adj_map.net');
    delete('adj_map.vec');
    delete('adj.clu');
    delete('adj.map');
    delete('adj.tree');

end
