function [A, B] = addDummyNodes(pieceA, pieceB)

    nPieceA = length(pieceA);
    nPieceB = length(pieceB);
    
    % If the number of vertices from each graph are not equal in the
    % cluster, add dummy vertices for the smaller graph
    if nPieceA < nPieceB
        % To handle matching dummy vertices, we turn true nonedges in the
        % adjacency matrices into -1's
        A = 2*pieceA - ones(nPieceA);
        B = 2*pieceB - ones(nPieceB);
        diff = nPieceB - nPieceA;
        Z12 = zeros(nPieceA, diff);
        Z22 = zeros(diff);
        A = [A, Z12; Z12', Z22]; % dummy vertices are disconnected
    elseif nPieceA > nPieceB
        % To handle matching dummy vertices, we turn true nonedges in the
        % adjacency matrices into -1's
        A = 2*pieceA - ones(nPieceA);
        B = 2*pieceB - ones(nPieceB);
        diff = nPieceA - nPieceB;
        Z12 = zeros(nPieceB, diff);
        Z22 = zeros(diff);
        B = [B, Z12; Z12', Z22]; % dummy vertices are disconnected
    else
        A = pieceA;
        B = pieceB;
    end
end