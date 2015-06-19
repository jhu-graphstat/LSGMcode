
%% A function for Dan S to test that things kind of work

try
    rng(job);
    fprintf('Job is #%i\n',job);
catch
    fprintf('No job specified. Defaulting to job is #1\n');
    rng(40);
end

%% Start in the main git repo directory
clear

r = rng;
job = double(r.Seed);

addpath(genpath('.'))

% %% Load the data
% 
% % YOU NEED TO CHANGE THIS IF YOU ARE NOT DAN SUSSMAN
% 
% dir_04 = '~/Dropbox/Data/twitter_xdata_2015-05-16/2014_04/';
% node_fn_04 = [dir_04 'twitter_prune_w50000.nodes.txt'];
% edge_fn_04 = [dir_04 'twitter_prune_w50000.edges.txt'];
% [node_04, adj_04] = readTwitter(node_fn_04,edge_fn_04);
% 
% dir_05 = '~/Dropbox/Data/twitter_xdata_2015-05-16/2014_05/';
% node_fn_05 = [dir_05 'twitter_prune_w50000.nodes.txt'];
% edge_fn_05 = [dir_05 'twitter_prune_w50000.edges.txt'];
% [node_05, adj_05] = readTwitter(node_fn_05,edge_fn_05);
% 
% %% Get common vertices
% [node_common, adj_04_common, adj_05_common] = ...
%     commonVerticesTwitter(node_04,adj_04, node_05, adj_05);
% 
% %% Reduce to only users
% node_common.is_hash = cellfun(@(name) name(1)=='#', node_common.name);
% 
% adj_04_user = adj_04_common(~node_common.is_hash,~node_common.is_hash);
% adj_05_user = adj_05_common(~node_common.is_hash,~node_common.is_hash);

%% Load the data fast

load data/twitterSmall
fprintf('Loaded Data\n');
%% Match common user graphs

nIter = 20;
start = 'convex';
nSeeds = ceil(job/10)*50+50;

normalize = @(x) x/max(x);

weight_table = {...
    'id',           @(x)x; ...
    'rank',         'rank';...
    'log(x)+1',     @(x)log(x)+1;...
    'sqrt',         @sqrt;...
%     'id norm',      normalize;...
%     'log norm',     {@(x)log(x)+1,normalize};...
%     'sqrt norm',    {@sqrt,normalize};...
    'binary',       @(x)double(x>0)};

[nWeights,~] = size(weight_table);
    

error = ones(nWeights,nIter);
t = zeros(nWeights,nIter);
weight = cell(nWeights,nIter);
    
totalRuns = nIter*nWeights;
jobId = double(job)*ones(totalRuns,1);
numSeeds = double(nSeeds)*ones(totalRuns,1);
%%
for i=1:nIter
    fprintf('\nIter #%i\n',i);
    %% Set Seeds
    [nUser,~] = size(adj_04_user);
    seeds = sort(datasample(1:nUser,nSeeds, 'Replace',false));
    nonSeeds = true(1,nUser);
    nonSeeds(seeds) = false;
    nonSeedsIdx = find(nonSeeds);
    
    %% Loop over weights and match
    for j=1:nWeights
        weightName = weight_table{j,1};
        weightFunc = weight_table{j,2};
        fprintf(weightName)
        tic;
        [coor,~] = seedgraphmatchell2(reweight(adj_04_user,weightFunc),...
                                      reweight(adj_05_user,weightFunc),...
                                      seeds, start);
        t(j,i) = toc;    
        error(j,i) = 1-mean(coor(nonSeeds)==nonSeedsIdx);
        weight{j,i} = weightName;
        fprintf('Error=%f. Time=%f \n',error(j,i),t(j,i));
        save(sprintf('results/weight_test_v02_corr_job%i_mc%i_weight_%s',...
            job, i, weightName), 'seeds','nonSeeds', 'coor');
    end

    
    %% Store in a table ...
    weight_col= reshape(weight,[nIter*nWeights,1]);
    error_col = reshape(error,[nIter*nWeights,1]);
    t_col = reshape(t,[nIter*nWeights,1]);

    results = table(weight_col, error_col, t_col,numSeeds,jobId);
    results.Properties.VariableNames = {'weight' 'error' 'time' 'numSeeds','jobId'};
    save(sprintf('results/weight_test_v02_job%i',job),'results');
end
 

%% save results
if ~exist('job')
    allRes = table();
    for job=1:40
        load(sprintf('results/weight_test_v02_job%i',job));
        allRes = [allRes; results];
    end
    %%
    allRes = allRes(~cellfun(@isempty, allRes.weight),:);
    weightingTestSummary = grpstats(allRes,{'weight','numSeeds'},{'mean'},...
                     'DataVars',{'error','time'})
    writetable(allRes);
end
