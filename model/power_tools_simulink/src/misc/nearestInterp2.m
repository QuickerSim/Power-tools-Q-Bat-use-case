function target_values = nearestInterp2(x, y, v, x_target, y_target)
coder.allowpcode('plain')
x_half_expnaded = [x(1, :); x; x(end, :)];
x_expanded = [-1e2 * ones(size(x_half_expnaded, 1), 1), ...
    x_half_expnaded, ...
    1e2 * ones(size(x_half_expnaded, 1), 1)];

y_half_expnaded = [y(:, 1), y, y(:, end)];
y_expanded = [-1e2 * ones(1, size(y_half_expnaded, 2)); y_half_expnaded; 1e2 * ones(1, size(y_half_expnaded, 2))];

v_expanded = nan(size(v) + [2, 2]);
v_expanded( 2 : end - 1, 2 : end - 1 ) = v;
v_expanded(1,  2 : end - 1) = v(1, :);
v_expanded(end,  2 : end - 1) = v(end, :);
v_expanded(2 : end - 1,  1) = v(:, 1);
v_expanded(2 : end - 1,  end) = v(:, end);
v_expanded(1, 1) = v(1, 1);
v_expanded(1, end) = v(1, end);
v_expanded(end, end) = v(end, end);
v_expanded(end, 1) = v(end, 1);


target_values = interp2( ...
    x_expanded, y_expanded, v_expanded, x_target, y_target, 'linear');
end