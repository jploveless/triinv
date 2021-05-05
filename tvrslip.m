function x = tvrslip(G,W,D,Patches,lambda,varargin);

%TVRslip(G,W,D,Patches,lambda,lock,nneg) 
%Lock subjects the inversion to locking constraints on the up- and/or
%downdip extent of the triangular mesh. Lock is a 1x3 logical in the form
%[updip downdip lateral] that chooses which edges to lock.

% Set up weighted matrices
A = sparse((W^(1/2))*G); 
b = sparse((W^(1/2))*D);

% Indices of distinct faults
nTri = sum(Patches.nEl);
ends = cumsum(Patches.nEl);
begs = [1; ends(1:end-1) + 1];

% Parse optional input arguments
if nargin < 6 % If no optional arguments were specified, 
   lock = false(length(Patches.nEl), 3); % No edge constraints
   nneg = zeros(length(Patches.nEl), 2); % No sign constraints
elseif nargin > 5 % If optional arguments were specified,
   for i = 1:length(varargin)
      if size(varargin{i}, 2) == 3 % Edge locking is 3 columns
         lock = double(varargin{i});
      elseif size(varargin{i}, 2) == 2 % 2-column array is sign constraint
         nneg = varargin{i};
      elseif size(varargin{i}, 2) == nTri % Big matrix is diff matrix
         Diff = varargin{i};
      end
   end
end   
%
%Apply lambdas across elements
%Lambda = zeros(2*nTri, 1);
%if numel(lambda) == 1
%   Lambda = Lambda + lambda;
%elseif numel(lambda) == numel(Patches.nEl)   
%   for i = 1:length(lambda)
%      Lambda(2*begs(i)-1:2*ends(i)) = lambda(i);
%   end
%else 
%   error('lambda should be a scalar applied to all faults or a vector with one value per fault.')   
%end
  
% Handle edge locking

% Blank index arrays
toplock = zeros(nTri, 2);
botlock = toplock;
sidlock = toplock;

% Loop through each fault
if size(lock, 1) == 1
   lock = repmat(lock, length(Patches.nEl), 1);
end
for i = 1:length(Patches.nEl)
   if sum(lock(i, :)) ~= 0
      edges = edgeelements(Patches.c, Patches.v(begs(i):ends(i), :));
      toplock(begs(i):ends(i), :) = lock(i, 1).*repmat(edges.top, 1, 2);
      botlock(begs(i):ends(i), :) = lock(i, 2).*repmat(edges.bot, 1, 2);
      sidlock(begs(i):ends(i), :) = lock(i, 3).*repmat((edges.s1|edges.s2), 1, 2);
   end      
end

toplock = logical(stack2(toplock));
botlock = logical(stack2(botlock));
sidlock = logical(stack2(sidlock));

% Handle sign constraints
spos = false(nTri, 2);
sneg = spos;
dpos = spos;
dneg = spos;

if exist('nneg', 'var')
   for i = 1:size(nneg, 1)
      if nneg(i, 1) == 1
         spos(begs(i):ends(i), 1) = true;
      elseif nneg(i, 1) == -1   
         sneg(begs(i):ends(i), 1) = true;
      end
      if nneg(i, 2) == 1
         dpos(begs(i):ends(i), 2) = true;
      elseif nneg(i, 2) == -1   
         dneg(begs(i):ends(i), 2) = true;
      end
   end
end

spos = stack2(spos);
sneg = stack2(sneg);
dpos = stack2(dpos);
dneg = stack2(dneg);
% Create diff matrix if not passed as input
if ~exist('Diff', 'var')
   Diff = MakeDiffMatrix_mesh2d(Patches);
end

% Estimate slip
n = size(A, 2);
cvx_begin quiet
variable x(n)
% Minimization statement
minimize( norm(A*x-b,2) + (lambda).*norm(Diff*x,1) )
subject to 
    % Sign constraints
    x(dpos) >= 0;
    x(dneg) <= 0;
    x(spos) >= 0;
    x(sneg) <= 0;
    % Edge constraints
    x(toplock) == 0;
    x(botlock) == 0;
    x(sidlock) == 0;
cvx_end

end
