function P = topKMatching(A, B, seeds, alpha, numRestarts, center, similarity, similarityWeight)

[A_full, B_full] = addDummyNodes(A, B);

roundSGM = @(A, B, s, topK, init) sgmSimilarity(A, B, s, topK, init, similarity, similarityWeight, true);

P = alphaSpokeGraphMatching(A_full, B_full, seeds, true, alpha, numRestarts, center, roundSGM);