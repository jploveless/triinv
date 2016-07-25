function [u, g, w, We, d] = triinvx(s, p, beta, k, verttrim, nneg, varargin)
%
% TRIINVX   Inverts surface displacements for slip on triangular mesh in Cartesian space.
%    TRIINV(S, T, BETA, KERNEL, VT) inverts the displacements/velocities contained 
%    in the .sta.data file S for slip on the triangular dislocation mesh defined
%    in the file T, which can be in .mat or .msh format (see below), subject to the
%    Laplacian smoothing constraint whose strength is defined as BETA.  The KERNEL 
%    argument allows for specification of an existing matrix containing the triangular
%    dislocation element partial derivatives, or the name of a file to which the partials
%    will be saved. VT is a flag indicating whether vertical velocities should be ignored
%    (VT = 1) or considered (VT = 0). 
%
%    TRIINV(S, T, BETA, KERNEL, LOCK) also subjects the inversion to locking constraints
%    on the up and/or downdip extents of the triangular mesh.  LOCK is a vector
%    containing N or 3N elements, where N is the number of distinct patches in T.
%    For example, to lock the updip row of elements while placing no constraints
%    on the downdip extent, or lateral elements, LOCK = [1 0 0].  It is important
%    to note that the locking procedure finds the up- and downdip elements assuming
%    that the fault has been constructed using a depth contour constraint.  That is,
%    updip elements are defined as those having 2 nodes at the maximum z-value (i.e.,
%    shallowest), while downdip elements have 2 nodes at the minimum (most negative)
%    z-value.  Lateral elements are simply edge-lining elements that belong to neither
%    the up- nor downdip set.
%
%    *** Important notes about the KERNEL argument: ***
%    If you specify a string for KERNEL, the program will check for its existence, and 
%    load it if it exists. If not, the elastic partials will be calculated and saved to
%    KERNEL. No internal checks are made to assure that an existing KERNEL corresponds 
%    to the present station-source geometry, except that the inversion will fail if the 
%    loaded kernel does not have the same dimensions as the present configuration 
%    of S and T.  To assure calculation of a "fresh" kernel, either delete the
%    existing kernel or provide a unique file name.
%
%    *** Notes about input argument T: ***
%    T can either be the path to a .msh file as written by the open-source meshing 
%    program Gmsh, or a .mat file containing two variables "c" and "v".  c is an n-by-3 
%    array whose columns contain the X, Y, Z coordinates of the n nodes in a triangular
%    mesh.  Z coordinates should be negative, i.e. a depth of 15 km would be given as -15.
%    v is an m-by-3 array containing the indexes of the nodes that comprise each triangle.
%    For example, if element 1 of the mesh is made up of nodes 10, 17, and 103, the first
%    line of V should be [10, 17, 103].  
%
%    U = TRIINV(...) returns the estimated slip (rate) to vector U.  U is structured
%    so that strike and dip slip estimates are stacked, i.e. 
%       [Us1;
%        Ud1;
%        Us2;
%        Ud2;
%    U =  . 
%         .
%         .
%        Usm;
%        Udm];

% Check station coordinates
%[s.x, s.y, s.z] = deal(s.lon, s.lat, 0*s.lat);
if ~isfield(s, 'z')
   s.z                              = 0*s.x;
end

% Add vertical velocity and uncertainty, if need be
if ~isfield(s, 'upVel')
   s.upVel                          = 0*s.eastVel;
end
if ~isfield(s, 'upSig')
   s.upSig                          = ones(length(s.lon), 1);
end

if ischar(p) % if the triangular mesh was specified as a file...
   % load the patch file...
   p                                = ReadPatches(p);
else
   % Check to see that the necessary fields were provided
   if sum(isfield(p, {'c', 'v', 'nEl', 'nc'})) ~= 4
      error('Missing one of more required fields in patch structure.', 'missingPatchf');
   end
end

% Process patch coordinates
p                                   = PatchCoordsx(p);
% Identify dipping and vertical elements
tz                                  = 3*ones(p.nEl, 1); % By default, all are vertical (will zero out second component of slip)
tz(abs(90-p.dip) > 1)               = 2; % Dipping elements are changed

% Check the length of beta and repeat if necessary
if numel(beta) ~= numel(p.nEl)
   beta                             = repmat(beta, numel(p.nEl), 1);
end

% Check for existing kernel
if ischar(k)
   if exist(k, 'file')
      load(k)
   else
      % Calculate the triangular partials
      g                             = GetTriCombinedPartialsx(p, s, [1 0]);
      save(k, 'g');
   end
else
   g                                = k;
end

% Adjust triangular partials
triS                                = [1:length(tz)]'; % All include strike
triD                                = find(tz(:) == 2); % Find those with dip-slip
triT                                = find(tz(:) == 3); % Find those with tensile slip
colkeep                             = setdiff(1:size(g, 2), [3*triD-0; 3*triT-1]);
gt                                  = g(:, colkeep); % eliminate the partials that equal zero

% Trim rows vertical corresponding to vertical displacements?
if verttrim == 1
   rowkeep                          = sort([[1:3:3*numel(s.x)], [2:3:3*numel(s.x)]]');
else
   rowkeep                          = 1:3*numel(s.lon);
end   

gt                                  = gt(rowkeep, :);

% Lock edges?
if nargin > 6
   Command.triEdge                  = varargin{:};
else
   Command.triEdge                  = 0;
end

if sum(Command.triEdge) ~= 0
   Ztri                             = ZeroTriEdges(p, Command);
else
   Ztri                             = zeros(0, size(gt, 2));
end

% Spatial smoothing
if sum(beta)  
   [w, be]                          = trismoothmp(p, s, gt, beta);
else
   w                                = zeros(0, 3*size(p.v, 1));
end

% Weighting matrix
we                                  = stack3(1./[s.eastSig, s.northSig, s.upSig].^2);
we                                  = we(rowkeep);
we                                  = [we ; be.*ones(size(w, 1), 1)]; % add the triangular smoothing vector
we                                  = [we ; 1e10*ones(size(Ztri, 1), 1)]; % add the zero edge vector
We                                  = diag(we); % assemble into a matrix

% Assemble the Jacobian...
G                                   = [gt; w; Ztri];
% ...and the data vector
d                                   = stack3([s.eastVel, s.northVel, s.upVel]);
d                                   = d(rowkeep);
d                                   = [d; zeros(size(w, 1), 1); zeros(size(Ztri, 1), 1)];

  dcov = stack3([s.eastSig, s.northSig, s.upSig]).^2; % Data covariance
  dcov = diag(dcov(rowkeep));
  dcho = chol(dcov); % Cholesky decomposition of data covariance
  dcInv = inv(dcho');
  Data = stack3([s.eastVel, s.northVel, s.upVel]);
  Data = Data(rowkeep);
  dataW = dcInv*Data; % Weighed data
  gW = dcInv*gt; % Weighted partials
  % Assemble matrices 
  bigD = [dataW; zeros(size(w, 1) + size(Ztri, 1), 1)]; % Data vector, augmented with smoothing and edge constraints
  bigG = [gW; diag(be)*w; 1e10.*Ztri]; % Partials, augmented with smoothing and edge partials


if nneg == 0
% Backslash inversion
	% Carry out the inversion
%	u                                   = (G'*We*G)\(G'*We*d);
   u                                   = bigG\bigD;
else
% Use non-negative solver

   % Do the inversion

   options = optimoptions('lsqlin', 'tolfun', 1e-25, 'maxiter', 1e5, 'tolpcg', 1e-3, 'PrecondBandWidth', Inf);
   u = lsqlin(bigG, bigD, [], [], [], [], repmat([-.1; 0], size(bigG, 2)/2, 1), repmat([.1; Inf], size(bigG, 2)/2, 1), [], options);
%   u = lsqlin(G'*We*G, G'*We*d, [], [], [], [], repmat([-.1; 0], size(G, 2)/2, 1), repmat([.1; Inf], size(G, 2)/2, 1), [], options);
end


% Pad the estimated slip with zeros, corresponding to slip components that were not estimated
U                                   = zeros(3*numel(p.xc), 1);
U(colkeep)                          = u;
u                                   = U;
