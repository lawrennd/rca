function visualiseNetwork(S, edgeColor, nodeLabels)

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
arcs = linspace(pi/2, pi/2 - 2*pi, nNodes+1)';
nodePositions = [cos(arcs(1:end-1)) sin(arcs(1:end-1))];
plot(nodePositions(:,1), nodePositions(:,2),'.k', 'MarkerSize',1)
ylim([-1.5 1.5])
xlim([-1.5 1.5])
daspect('manual')
% set(gca,'YTick',[], 'XTick',[])
axis image off
fontsize = 18;

if nargin < 2
    edgeColor = 'k';
end
for i = 1:nNodes
    for j = 1:nNodes
        if S(i,j) ~= 0
            hold on, fill(nodePositions([i j],1), nodePositions([i j],2), edgeColor, 'linewidth', 2);
        end
    end
    arcs = linspace(0, 2*pi, 360)';
    nChars = numel(nodeLabels{i});
    nodeBubble = [nChars*.05*cos(arcs(1:end-1)) .05*sin(arcs(1:end-1))];
%     hold on, plot(nodePositions(i,1)+nodeBubble(:,1), nodePositions(i,2)+nodeBubble(:,2),'k')
    hold on, fill(nodePositions(i,1)+nodeBubble(:,1), nodePositions(i,2)+nodeBubble(:,2), [1 1 1], 'EdgeColor', 'k'),
end
htext = text(nodePositions(:,1)-.1, nodePositions(:,2), nodeLabels);
set(htext, 'fontsize', fontsize, 'fontweight', 'bold', 'color', 'k');

