function target_values = nearestInterp1(argument, values, target_arguments)
coder.allowpcode('plain')
target_values = interp1( ...
    [-1e12; argument; 1e12], [values(1); values; values(end)], target_arguments, 'linear');
end