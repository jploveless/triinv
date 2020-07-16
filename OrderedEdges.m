function elo = OrderedEdges(c, v);
% OrderedEdges returns the edges 
el = boundedges(c, v);
elo = [el(1, :)];

for i = 2:length(el);
	[next, col] = find(el == elo(i-1, 2)); % find all of the boundary lines containing the second entry of the current ordered boundary line
	if numel(next) > 2
   	   [~, n] = setdiff(sum(el(next, :), 2), sum(elo, 2), 'rows'); % choose that which is not the current boundary line
    else
       n = find(sum(el(next, :), 2) ~= sum(elo(i-1, :), 2)); 
    end
    if isempty(n) % Handle holes in the mesh
       keyboard
    end
	next = next(n(1)); col = col(n(1)); 
	elo = [elo; [el(next, col) el(next, setdiff([1 2], col))]]; % order the endpoints of the next boundary line
end
elo = elo';