function h = plotDataSet(h, X, y, t)
colors = {'r', 'k', 'b'};
% styles = {'o', 's', '^'};

for j = 1:3
    scatter(h, X(y == j, 1), X(y == j, 2), 28, 'o', colors{j}, 'filled');
    hold(h, 'on');
end
axis(h, 'square');
title(h, t, 'FontSize', 10);