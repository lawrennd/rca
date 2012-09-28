function visualiseNetwork(S, varargin)

% VISUALISENETWORK Visualise network graph.
%
% FORMAT
% DESC Visualises the network represented by the connectivity matrix S.
% An edge is drawn between nodes i and j, if S(i,j) is non-zero.
%
% ARG S : the connectivity matrix.
%
% The function also accepts a variable length input argument list, in the
% format ['arg1', arg1 , 'arg2', arg2, ...].
% Available options are:
%
%   'edgeColor' : the string of the color used for plotting the edges.
%   If edgeColor is a matrix of the same size as S, then every edge is
%   colored based on edgeColor(i,j) and the current colormap.
%
%   'nodeLabels' : the cell array of labels to annotate the nodes.
%
%   'nodePositions' : the size(S,1) x 2 matrix of the node coordinates.
%   Note that for optimum effect, the coordinates are automatically
%   normalised within the [-1 1] interval, unless the option normalise is
%   set to false.
%
%   'curv' : the edge curvature. Zero for straight edges.
%
%   'normalise' : the boolean indicator for normalising the coordinates of
%   the nodes positions.
%
%   'fontSize' : the node label font size.
%
%   'nodeSize' : the node size as a fraction of the label size. If nodeSize
%   is 0 then no nodes are rendered.
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

nNodes = size(S,1);

arcs = linspace(pi/2, pi/2 - 2*pi, nNodes+1)';                              % Default node circle formation.
nodePositions = [cos(arcs(1:end-1)) sin(arcs(1:end-1))];
edgeColor = repmat({'k'},size(S));                                          % Default edge color.
nodeLabels = num2cell(1:nNodes);                                            % Default node labels.
fs = 12;                                                                    % Default font size.
ns = 1;                                                                     % Default node size.

for k = 1:2:length(varargin)
    switch varargin{k}
        case 'nodePositions'; nodePositions = varargin{k+1};
            [tf,loc] = ismember('normalise', varargin(1:2:length(varargin)));
            if tf
                normalise = varargin{2*loc};
                if normalise
                    nodePositions = nodePositions ./ repmat(max(nodePositions), nNodes, 1); % Normalise coordinates.
                end
            end
        case 'curv'; curv = varargin{k+1};
        case 'nodeLabels'; nodeLabels = varargin{k+1};
        case 'edgeColor';
            if numel(varargin{k+1}) == numel(S)
                edgeColor = varargin{k+1} - diag(diag(varargin{k+1}));
                if isnumeric(varargin{k+1})
                    edgeColorCarray = num2cell(edgeColor);
                    colorRange = unique(edgeColor);
                    cmap = jet(length(colorRange));
                    for i = 1:nNodes
                        for j = 1:nNodes
                            if i < j && S(i,j) ~= 0
                                edgeColorCarray{i,j} = cmap(colorRange == edgeColor(i,j), :);
                            end
                        end
                    end
                else iscell(varargin{k+1})
                    edgeColorCarray = varargin{k+1};
                end
            elseif ischar(varargin{k+1}) || (isnumeric(varargin{k+1}) && length(varargin{k+1})==3)
                edgeColorCarray = repmat(varargin(k+1), size(S));
            end
        case 'fontSize'; fs = varargin{k+1};
        case 'nodeSize'; ns = varargin{k+1};
    end
end

plot(nodePositions(:,1), nodePositions(:,2),'.k', 'MarkerSize',.1), hold on
if normalise
    ylim([-.1 1.1]), xlim([-.1 1.1])
end

for i = 1:nNodes
    for j = 1:nNodes
        if i < j && S(i,j) ~= 0
            % fill(nodePositions([i j],1), nodePositions([i j],2), 'k', 'edgecolor',edgeColorCarray{i,j}, 'linewidth',1);
            % line(nodePositions([i j],1), nodePositions([i j],2), 'LineStyle','-', 'color',edgeColorCarray{i,j}, 'linewidth',.5);
            [xc, yc] = curve(nodePositions([i j],1), nodePositions([i j],2), curv);
            plot(xc,yc,'LineStyle','-', 'color', edgeColorCarray{i,j}, 'linewidth',.5)
        end
    end
end

arcs = linspace(0, 2*pi, 360)';
% nChars = numel(nodeLabels{i});
if ~isempty(nodeLabels)
    htext = text(nodePositions(:,1), nodePositions(:,2), nodeLabels);           % Plot node labels.
    set(htext, 'fontsize', fs, 'fontweight', 'bold', 'color', 'k', 'HorizontalAlignment','center', 'verticalAlignment', 'middle');
end

% Plot nodes.
% nodeBubble = [nChars*(fs/2000)*cos(arcs(1:end-1)) (fs/400)*sin(arcs(1:end-1))];
if ns > 0
    for i = 1:nNodes
        txtSize = handle(htext(i)).Extent(3:4);
        nodeBubble = [.65*ns*txtSize(1)*cos(arcs(1:end-1)) .5*ns*txtSize(2)*sin(arcs(1:end-1))];
        fill(nodePositions(i,1)+nodeBubble(:,1), nodePositions(i,2)+nodeBubble(:,2), [1 1 1], 'EdgeColor', 'k')
    end

    allChildren = get(gca,'Children');
    textChildren = findobj(allChildren,'Type','text');
    set(gca,'Children',[textChildren; setdiff(allChildren,textChildren)]);
end
daspect('manual'), axis image off
h = colorbar; colormap jet
set(gca, 'CLim', [min(min(edgeColor)), max(max(edgeColor))])
x1 = get(h, 'Position');    x1(3) = x1(3)/5;
x2 = get(gca, 'Position');
set(h, 'Position', x1)
set(gca, 'Position', x2)
end

function [xc, yc] = curve(x, y, pct)                                        % Make a bent line between two points.
    if diff(x) == 0
        midx = x(1) + pct;  midy = sum(y)/2;
    elseif diff(y) == 0
        midx = sum(x)/2;    midy = y(1) + pct;
    else
        midx = sum(x)/2 + pct;  midy = sum(y)/2 + pct;
    end
    xc = spline([1 3 2],[x; midx], 1:1/10:3);
    yc = spline([1 3 2],[y; midy], 1:1/10:3);
end
