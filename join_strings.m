function log_vars = join_strings(logs)
    log_vars = '[';
    for i = 1:length(logs)
        if i > 1
            log_vars = [log_vars, ', ', logs{i}];
        else
            log_vars = [log_vars, logs{i}];
        end
    end
    log_vars = [log_vars, ']'];
end