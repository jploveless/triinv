% test_dip   Example estimates using triinvx

%### Basic estimate
% Load example inputs
load teststation-dip.mat
p = ReadPatches('testmesh-dip.msh');
load inslip-dip.mat
% Run models
u = triinvx(p, s, 100); % Simple estimation using beta = 100
[u, pred, gsave] = triinvx('testmesh-dip.msh', s, 100); % Same model but specifying the fault file as an input string and requesting predicted displacement field
% Plot synthetic slip that was used to generate the constraining displacements
meshview(p.c, p.v, utrue(2:3:end));
title('Input slip distribution');
caxis([-1 1]);
colormap(bluewhitered);
% Visualize estimated dip slip in a new figure
meshview(p.c, p.v, u(2:3:end));
title('Basic slip estimate, beta = 100');
caxis([-1 1]);
colormap(bluewhitered);
% Add observed and estimated displacements
hold on
scale_factor = 50; % Manual scaling factor
quiver(s.x, s.y, scale_factor*s.eastVel, scale_factor*s.northVel, 0); % Plot data as manually scaled vectors
quiver(s.x, s.y, scale_factor*pred(1:3:end), scale_factor*pred(2:3:end), 0); % Plot predictions as manually scaled vectors

%### Estimate with different smoothing weight, reusing the partial derivatives output in line 10
[u_500, pred_500] = triinvx(p, s, 500, 'partials', gsave);
meshview(p.c, p.v, u_500(2:3:end));
title('Basic slip estimate, beta = 500');
caxis([-1 1]);
colormap(bluewhitered);

%### Estimate enforcing reverse slip

[u_rev, pred_rev] = triinvx(p, s, 100, 'nneg', [0 1], 'partials', gsave);
meshview(p.c, p.v, u_rev(2:3:end));
title('Force reverse slip, beta = 100');
caxis([-1 1]);
colormap(bluewhitered);

%### Estimate enforcing reverse slip and no slip at fault edges

[u_rev_noedge, pred_rev_noedge] = triinvx(p, s, 100, 'nneg', [0 1], 'lock', [1 1 1], 'partials', gsave);
meshview(p.c, p.v, u_rev_noedge(2:3:end));
title('Force reverse slip, no edge slip, beta = 100');
caxis([-1 1]);
colormap(bluewhitered);

%### Estimate using total variation regularization, outputs requested as structures

[u_tvr, pred_tvr] = triinvx(p, s, 1e-2, 'tvr', true, 'outstruct', true, 'partials', gsave);
meshview(p.c, p.v, u_tvr.dip);
title('TVR, beta = 0.01');
caxis([-1 1]);
colormap(bluewhitered);

%### Estimate using total variation regularization, forcing reverse slip and no slip at edges, outputs requested as structures

[u_tvr_rev_noedge, pred_tvr_rev_noedge] = triinvx(p, s, 1e-2, 'tvr', true, 'lock', [1 1 1], 'nneg', [0 1], 'outstruct', true, 'partials', gsave);
meshview(p.c, p.v, u_tvr_rev_noedge.dip);
title('TVR, force reverse slip, no edge slip, beta = 0.01');
caxis([-1 1]);
colormap(bluewhitered);
