function [u, pred, gsave] = triinvx(p, s, beta, varargin)
%
% TRIINVX   Inverts surface displacements for slip on triangular mesh in Cartesian space.
%    TRIINV(P, S, BETA) inverts the displacements/velocities contained 
%    in the structure S for slip on the triangular dislocation mesh defined
%    in the file P, which can be in .mat or .msh format (see below), subject to the
%    Laplacian smoothing constraint whose strength is defined as BETA. 
%
%    The input structure S should contain coordinates as S.x, S.y, and S.z; 
%    displacements/velocities as S.eastVel, s.northVel, s.upVel; uncertainties as
%    s.eastSig, s.northSig, and s.upSig. If stresses are to be included in the 
%    inversion, the coordinates should be specified as S.xs, S.ys, and S.zs; tensor 
%    components at these coordinates should be fields S.sxx, S.syy, S.szz, S.sxy, 
%    S.sxz, S.syz, with tensor component uncertainties of S.sxxs, S.syys, S.szzs, 
%    S.sxys, S.sxzs, S.syzs. These fields can be generated from a 
%    6*nCoordinates-by-1 vector using the function makestressfields.m. 
%
%    Optional input arguments can be specified using 'parameter', value pairs:
%
%    TRIINV(..., 'partials', G) enables specification of a pre-calculated matrix, 
%    G, of partial derivatives relating slip on triangular dislocation elements to 
%    displacement at observation coordinates. G should be untrimmed in both the 
%    rows in columns.
%
%    TRIINV(..., 'lock', EDGES) subjects the inversion to locking (zero slip) 
%    constraints on the updip, downdip, and/or lateral edges of the triangular mesh.  
%    EDGES is a 3-column array with non-zero values corresponding to the edges(s) 
%    to be locked. For example, EDGES = [1 0 1] imposes no-slip conditions on the updip
%    and lateral edges of the fault (columns 1 and 3, respectively), but not on the 
%    downdip edge (column 2). For structures T that contain multiple distinct faults
%    (i.e., numel(T.nEl) > 1), specify EDGES as a single row to apply the constraints to 
%    all faults, or include a unique row for each fault. 
%
%    TRIINV(..., 'dcomp', DCOMPONENTS) allows specification of which components of the 
%    displacement (velocity) observations to use. DCOMPONENTS is a 3-element vector with
%    non-zero entries corresponding to the components to consider as observations in the 
%    inversion. For example, to use only vertical displacement data, specify:
%    DCOMPONENTS = [0 0 1]. The default behavior is DCOMPONENTS = [1 1 1] so that the 
%    X, Y, and Z components are all considered as constraining data in the inversion. 
% 
%    TRIINV(..., 'scomp', SCOMPONENTS) allows specification of which components of the 
%    stress (stress rate) observations to use. SCOMPONENTS is a 6-element vector with
%    non-zero entries corresponding to the components to consider as observations in the 
%    inversion. For example, to use only shear stress data, specify:
%    SCOMPONENTS = [0 0 0 1 1 1]. The default behavior is COMPONENTS = [1 1 1 1 1 1] 
%    so that all tensor components are considered as constraining data in the inversion. 
%
%    TRIINV(..., 'stresstype', devflag) allows specification of whether full or deviatoric
%    stress tensor components are given in the structure S. If devflag == 0 (or false), 
%    the stress tensor components are assumed to be full, and if devflag ~= 0 (or true), 
%    deviatoric components are assumed, and the Green's functions are modified by 
%    subtracting the mean stress Green's functions from the normal stress components. 
%
%    TRIINV(..., 'lame', MODULI) allows specification of Lame parameters defining the 
%    elasticity of the half-space. MODULI is a 2-element vector giving the Lame parameters
%    lambda and mu.
%
%    TRIINV(..., 'nneg', NONNEG) allows specification of which components of the 
%    slip distribution should be sign-constrained. NONNEG should be either a 2-element
%    vector that applies uniform constraints to all faults, or a numel(P.nEl)-by-2 array
%    that allows specification of different constraints for each fault. Each row of NONNEG
%    specifies constraints on the sign of [strike, dip] slip, with a 0 to place no 
%    constraint on the sign of slip, 1 to constrain slip in that direction to be positive,
%    and -1 to constrain slip to be negative. For example, to constrain dip-slip only to 
%    be positive, specify NONNEG = [0 1]. To constrain strike-slip to be negative and dip
%    slip to be positive, specify NONNEG = [-1 1]; 
%
%    TRIINV(..., 'tvr', TVRFLAG), where TVRFLAG = true, estimates slip using total
%    variation regularization rather than Laplacian smoothing. In this case, the 
%    regularization strength is still controlled by input argument BETA.
%
%    *** Notes about input argument P: ***
%    P can either be the path to a .msh file as written by the open-source meshing 
%    program Gmsh, or a .mat file containing two variables "c" and "v".  c is an n-by-3 
%    array whose columns contain the X, Y, Z coordinates of the n nodes in a triangular
%    mesh.  Z coordinates should be negative, i.e. a depth of 15 km would be given as -15.
%    v is an m-by-3 array containing the indexes of the nodes that comprise each triangle.
%    For example, if element 1 of the mesh is made up of nodes 10, 17, and 103, the first
%    line of V should be [10, 17, 103] (as returned from Matlab's delaunay.m routine).  
%
%    U = TRIINV(...) returns the estimated slip (rate) to vector U.  U is structured
%    so that strike, dip, and tensile slip estimates are stacked, i.e. 
%         U = [Us1, Ud1, Ut1, Us2, Ud2, Ut2, ..., Usm, Udm, Utm]';

% Parse optional inputs 
if nargin > 3
   for i = 1:2:numel(varargin)
      if startsWith(varargin{i}, 'partials')
         g = varargin{i+1};
      elseif startsWith(varargin{i}, 'lock')
         Command.triEdge = varargin{i+1};
      elseif startsWith(varargin{i}, 'dcomp')
         dcomp = varargin{i+1};
      elseif startsWith(varargin{i}, 'scomp')
         scomp = varargin{i+1};
      elseif startsWith(varargin{i}, 'stresstype')
         stresstype = logical(varargin{i+1});
      elseif startsWith(varargin{i}, 'lame')
         lame = varargin{i+1};
      elseif startsWith(varargin{i}, 'nneg')
         nneg = varargin{i+1};
      elseif startsWith(varargin{i}, 'tvr')
         tvr = varargin{i+1};
      elseif startsWith(varargin{i}, 'outstruct')
         outstruct = varargin{i+1};
      end
   end
end

% Check for existence of specified Lame parameters
if ~exist('lame', 'var')
   mu = 3e10; lambda = 3e10; % Defaults
else
   mu = lame(2); lambda = lame(1);
end

% (Re)calculate Poisson's ratio from Lame parameters
pr = lambda./(2*(lambda + mu));

% Check whether TVR should be used
if ~exist('tvr', 'var')
   tvr = false;
end

% Check station structure to make sure z coordinates exist; assume surface if not specified
if ~isfield(s, 'z')
   s.z                              = 0*s.eastVel;
end

% Station count
numsta = numel(s.eastVel);
if isfield(s, 'xs')
   numstress = numel(s.xs);
else
   numstress = 0;
end

% Check for specification of displacement components to use
if numsta > 0 && ~exist('dcomp', 'var')
   dcomp = [1 1 1]; % Default is to use all components
end

% Add template vertical velocity and uncertainty, if need be
if ~isfield(s, 'upVel')
   s.upVel                          = 0*s.eastVel;
   dcomp(3)                         = 0; % Make sure not to use these for inversion
end
if ~isfield(s, 'upSig')
   s.upSig                          = ones(numsta, 1);
end

% Check whether data coordinates are geographic or Cartesian
isgeog = false;
if isfield(s, 'lon')
   isgeog = true;
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
if isgeog
   p                                = PatchCoords(p);
else
   p                                = PatchCoordsx(p);   
end

% Identify dipping and vertical elements
tz                                  = 3*ones(sum(p.nEl), 1); % By default, all are vertical (will zero out second component of slip)
tz(abs(90-p.dip) > 1)               = 2; % Dipping elements are changed

% Check the length of beta and repeat if necessary
if numel(beta) ~= numel(p.nEl)
   beta                             = repmat(beta, length(p.nEl), 1);
end
% Apply beta(s) to every element
Beta                                = zeros(2*size(p.v, 1), 1); % Blank vector, 2 values for every element
ends                                = cumsum(2*p.nEl(:));
begs                                = [1; ends(1:end-1)+1];
for i = 1:length(p.nEl)
   Beta(begs(i):ends(i))            = beta(i);
end

% Check for existing kernel
if ~exist('g', 'var')
   % Calculate the triangular partials
   if isfield(s, 'sxx') % If stress components exist in station structure
      % Create new structure that concatenates displacement and stress coordinates
      [ss.x, ss.y, ss.z]            = deal([s.x; s.xs], [s.y; s.ys], [s.z; s.zs]);
      pdisp                         = [true(numsta); false(numstress)]; % Logical index of displacement coordinates
      pstr                          = ~pdisp; % Logical index of stress coordinates
      [gd, gs]                      = GetTriCombinedPartialsx(p, ss, [pdisp pstr], pr); % Calculate both partials
      gs                            = StrainToStressComp(gs', mu, lambda)'; % Convert strain to stress partials
   else
      if isgeog
         gd                         = GetTriCombinedPartials(p, s, [1 0], pr);
      else
         gd                         = GetTriCombinedPartialsx(p, s, [1 0], pr);
      end
      gs                            = zeros(0, size(gd, 2)); % Blank stress partial matrix
   end
else % Parse input partials
   gd                               = g(1:3*numsta, :); % First rows are displacements
   gs                               = g(end - (6*numstress - 1):end, :); % Last rows are stresses
end

% Check whether stress partials should be adjusted for deviatoric stresses
if exist('stresstype', 'var')
   if stresstype
      % Take average of normal stress components
      gsmeannorm = mean([gs(1:6:end, :), gs(2:6:end, :), gs(3:6:end, :)], 2);
      gsmeannorm = repmat(gsmeannorm, 1, size(gs, 2));
      % Place into array same size as gs
      gsm = zeros(size(gs));
      gsm(1:6:end, :) = gsmeannorm;
      gsm(2:6:end, :) = gsmeannorm;
      gsm(3:6:end, :) = gsmeannorm;
      % Correct stress partials by subtracting mean stress
      gs = gs - gsm;
   end
end

% Determine which rows of partial derivative matrix to keep.

% These rows correspond to displacement and stress components that are being considered.

% dcomp is the optional variable that controls displacements, and it should be a 
% 3-element vector with non-zero entries for those components that should be used.
% 
% If not specified as an input argument, it is set by default to use all 3 components.
dcomp                               = repmat(dcomp(:), numsta, 1); % Replicate for all stations
rowkeep                             = find(dcomp ~= 0);
gsave                               = [gd; gs]; % Full partials for output   
gd                                  = gd(rowkeep, :); % Trimmed displacement components

% scomp is the optional variable that controls stresses, and it should be a 
% 6-element vector with non-zero entries for those components that should be used.
if exist('scomp', 'var')
   scomp                            = repmat(scomp(:), numstress, 1);
   srowkeep                         = find(scomp ~= 0);
else
   srowkeep                         = 1:6*numstress;
end
gs                                  = gs(srowkeep, :);

% Stack adjusted partials matrices; this is basis for estimation
g                                   = [gd; gs]; 

% Adjust triangular partials
triS                                = [1:length(tz)]'; % All include strike
triD                                = find(tz(:) == 2); % Find those with dip-slip
triT                                = find(tz(:) == 3); % Find those with tensile slip
colkeep                             = setdiff(1:size(g, 2), [3*triD-0; 3*triT-1]);
gt                                  = g(:, colkeep); % eliminate the columns for unused slip components
%tres                                = sum(gt.^2, 1)'; % Sum of triangular partials, for weighting regularization

% Lock edges?
if ~exist('Command', 'var') % Command structure only exists if triEdge was specified as input argument
   Command.triEdge                  = false(1, 3);
end

if sum(Command.triEdge) ~= 0
   Ztri                             = ZeroTriEdges(p, Command);
else
   Ztri                             = zeros(0, size(gt, 2));
end

% Spatial smoothing
if sum(beta)  
   % Find neighbors
   share = SideShare(p.v);

   % Make the smoothing matrix
   sm = MakeTriSmoothAlt(share);
   sm = sm(colkeep, :);
   sm = sm(:, colkeep);
   sm = full(sm);
else
   sm                               = zeros(0, 3*size(p.v, 1));
end

% Weighting matrix
wd                                  = stack3(1./[s.eastSig, s.northSig, s.upSig].^2);
wd                                  = wd(rowkeep);
if isfield(s, 'sxx') % If stress data exist
   ws                               = stack6(1./[s.sxxs, s.syys, s.szzs, s.sxys, s.sxzs, s.syzs].^2);
   ws                               = ws(srowkeep);
   w                                = [wd ; ws]; % Combined weighting
else
   w                                = wd;
end
wc                                  = [w; Beta.*ones(size(sm, 1), 1)]; % add the triangular smoothing vector
%wc                                  = [w; Beta./tres]; % add the triangular smoothing vector
wc                                  = [wc; 1e10*ones(size(Ztri, 1), 1)]; % add the zero edge vector
Wc                                  = diag(wc); % assemble into a matrix

% Assemble the Jacobian...
G                                   = [gt; sm; Ztri];
% ...and the data vector
dd                                  = stack3([s.eastVel, s.northVel, s.upVel]);
dd                                  = dd(rowkeep);
if isfield(s, 'sxx') % If stress data exist
   ds                               = stack6([s.sxx, s.syy, s.szz, s.sxy, s.sxz, s.syz]);
   ds                               = ds(srowkeep);
   d                                = [dd; ds];
else
   d                                = dd;
end
dc                                  = [d; zeros(size(sm, 1), 1); zeros(size(Ztri, 1), 1)];

% Check inversion type
if ~exist('nneg', 'var') % Check for slip sense constraints
   nneg = 0;
end

if sum(abs(nneg(:))) == 0 % If no constraints, 
   if ~tvr
      % Use backslash inversion with Laplacian smoothing
      u                             = (G'*Wc*G)\(G'*Wc*dc);
   else
      % Use TVR optimization
      u                             = tvrslip(gt, diag(w), d, p, beta, logical(Command.triEdge));
   end      
else % If there are constraints on the sign
   if ~tvr 
      % Use built-in constrained solver
      u                             = constrainedslip(G, Wc, dc, nneg);
   else
      % Or pass constraints to TVR routine
      u                             = tvrslip(gt, diag(w), d, p, beta, logical(Command.triEdge), nneg);
   end
end

% Pad the estimated slip with zeros, corresponding to slip components that were not estimated
U                                   = zeros(3*sum(p.nEl), 1);
U(colkeep)                          = u;
u                                   = U;

% Predict displacements
pred                                = gsave*u; % gsave is the untrimmed matrix of partials

% Format output arguments as structure, if requested
if exist('outstruct', 'var')
   if outstruct
      us.strike = u(1:3:end);
      us.dip = 0*us.strike;
      us.tens = 0*us.strike;
      us.dip(tz == 2) = u(3*find(tz == 2)-1);
      us.tens(tz == 3) = u(3*find(tz == 3)-0);
      u = us; % Output slips
      Pred.eastVel = pred(1:3:3*numsta);
      Pred.northVel = pred(2:3:3*numsta);
      Pred.upVel = pred(3:3:3*numsta);
      stressvec = pred(end - (6*numstress - 1):end);
      if ~isempty(stressvec)
         Pred = makestressfields(stressvec);
      end
      pred = Pred; % Output predictions
   end
end
      
