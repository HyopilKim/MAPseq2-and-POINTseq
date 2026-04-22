function wheremaxtar(data,targets,lim)

[numRows, numCols] = size(data);

% Prepare the data for plotting
[x, ~] = meshgrid(1:numCols, 1:numRows); % Generate x and y grid
x = x(:); % Flatten the x grid
y_values = data(:); % Flatten the matrix (column-major order)

% Plot each value as a dot
if exist('lim','var')
    figure;scatter(x, y_values);xticks(1:numCols);xticklabels(targets);ylim([0 lim]);findfigs;
else
    figure;scatter(x, y_values);xticks(1:numCols);xticklabels(targets);findfigs;
end
