function clusters = find_clusters(sig_elec, neighbours)

electrodeList = {neighbours(sig_elec).label};  % example electrodes
neighbors = {neighbours(sig_elec).neighblabel};

n = numel(electrodeList);
nameToIdx = containers.Map(electrodeList, 1:n);
A = zeros(n, n);
for i = 1:n
    for j = 1:numel(neighbors{i})
        if isKey(nameToIdx, neighbors{i}{j})
            ni = nameToIdx(neighbors{i}{j});
            A(i, ni) = 1;
            A(ni, i) = 1;  % ensure symmetry
        end
    end
end

G = graph(A);
components = conncomp(G);  % gives cluster indices per electrode

numComponents = max(components);
clusters = cell(1, numComponents);
for i = 1:numComponents
    clusters{i} = electrodeList(components == i);
end