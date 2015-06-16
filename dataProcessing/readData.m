nodeFn1 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_1_nodes.txt';
edgeFn1 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_1_edges_relabeled.txt';
%nodeFn2 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_2_nodes.txt';
%edgeFn2 = '/home/jason/Documents/XDATA2015/summer-school/data/instagram_2_edges_relabeled.txt';

twitterNodeFn = '/home/jason/Documents/XDATA2015/summer-school/data/twitter_1_nodes.txt';
twitterEdgeFn = '/home/jason/Documents/XDATA2015/summer-school/data/toLSGM_006_twitter_edges_october_relabeled.txt';

B = readTwitter2(edgeFn1);
%B = readTwitter2(edgeFn2);
A = readTwitter2(twitterEdgeFn);

