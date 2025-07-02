% flattens cell arrays
function y = flatten(x)
    y = {};
    for i=1:length(x)
        if iscell(x{i})
            z = flatten(x{i});
            for j=1:length(z)
                y = [y, z{j}];
            end
        else
            y = [y, x{i}];
        end
    end
end