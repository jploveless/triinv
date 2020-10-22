function u = constrainedslip(G, We, d, nneg)
% constrainedslip   Estimates slip using lsqlin
%   u = constrainedslip(G, We, d, nneg) estimates slip using Green's functions
%   G, data uncertainty matrix We, data vector d, and sign constraints nneg. 
%   G, We, and d are are constructed within triinvx, and nneg is either a 
%   2-element vector that applies uniform constraints to all faults, or a 
%   numel(P.nEl)-by-2 array that allows specification of different constraints
%   for each fault. Each row of nneg specifies constraints on the sign of 
%   [strike, dip] slip, with a 0 to place no constraint on the sign of slip, 
%   1 to constrain slip in that direction to be positive, and -1 to constrain
%   slip to be negative. For example, to constrain dip-slip only to be positive, 
%   (reverse) specify nneg = [0 1]. To constrain strike-slip to be negative 
%   (dextral) and dip slip to be positive (reverse), specify nneg = [-1 1];


options = optimoptions('lsqlin'); % Use default linear least squares options
% Set bounds on sign of slip components
slipbounds = [-Inf 0; -Inf Inf; 0 Inf];
if size(nneg, 1) == 1 % If a single bound constraint is specified
   slipboundss = slipbounds(nneg(1) + 2, :);
   slipboundsd = slipbounds(nneg(2) + 2, :);
   lbound = repmat([slipboundss(1); slipboundsd(1)], size(G, 2)/2, 1);
   ubound = repmat([slipboundss(2); slipboundsd(2)], size(G, 2)/2, 1);
elseif size(nneg, 1) == numel(p.nEl) % If bounds are specified on each fault
   % Define element index ranges
   ends = cumsum(p.nEl(:));
   begs = [1; ends(1:end-1)+1];
   % Allocate space for slip bounds
   lbound = zeros(size(G, 2)/2, 2);
   ubound = lbound;
   % Loop through each fault to apply constraints to its elements
   for i = 1:numel(p.nEl)
      lbound(begs(i):ends(i), 1) = slipbounds(nneg(i, 1) + 2, 1); % Strike-slip lower bound
      lbound(begs(i):ends(i), 2) = slipbounds(nneg(i, 2) + 2, 1); % Dip-slip lower bound
      ubound(begs(i):ends(i), 1) = slipbounds(nneg(i, 1) + 2, 2); % Strike-slip lower bound
      ubound(begs(i):ends(i), 2) = slipbounds(nneg(i, 2) + 2, 2); % Dip-slip lower bound
   end
   % Reshape to column vectors
   lbound = stack2(lbound);
   ubound = stack2(ubound);
else
   error('Slip sign constraints should be specified as a 1-by-2 vector, applying the same constraints to all faults, or a numel(p.nEl)-by-2 array specifying unique constraints for each fault.')
end
u = lsqlin(G'*We*G, G'*We*d, [], [], [], [], lbound, ubound, [], options);