function metrics = absK_metrics_trim(metrics)
f = fieldnames(metrics);
n = metrics.isave;
for i = 1:numel(f)
    name = f{i};
    if isnumeric(metrics.(name)) && numel(metrics.(name)) >= n && ~strcmp(name,'isave')
        metrics.(name) = metrics.(name)(1:n,:);
    end
end
end
