function visualiseNetwork(S, edgeColor)

% VISUALISENETWORK Visualise network graph.
%
% FORMAT
% DESC Visualises the network represented by the matrix S. An edge is drawn
% between nodes i and j, if S_ij is non-zero.
% ARG S : the connection matrix.
% ARG EDGECOLOR : color used for plotting the edges.
%
% SEEALSO :
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

nNodes = size(S,1);
arcs = linspace(0, 2*pi, nNodes+1)';
nodePositions = [cos(arcs(1:end-1)) sin(arcs(1:end-1))];
plot(nodePositions(:,1), nodePositions(:,2),'.k', 'MarkerSize',1)
ylim([-1 1])
xlim([-1 1])
daspect('manual')
set(gca,'YTick',[], 'XTick',[])
if nargin < 2
    edgeColor = 'k';
end
for i = 1:nNodes
    for j = 1:nNodes
        if S(i,j) ~= 0
            line(nodePositions([i j],1), nodePositions([i j],2), 'LineWidth', 2, 'Color', edgeColor);
        end
    end
end
hold on
plot(nodePositions(:,1), nodePositions(:,2),'xk', 'MarkerSize',5)
hold off