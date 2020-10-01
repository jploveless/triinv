function w = MakeTriSmoothSimple(share)
%
% MakeTriSmoothSimple produces a scale-independent smoothing matrix.
%
% Inputs:
%   share		= n x 3 array of indices of the up to 3 elements sharing a side 
%					  with each of the n elements, from SideShare.m
%
% Outputs:
%   w		    = n x n sparse smoothing matrix
%

% 

n = size(share, 1);

% allocate space for the smoothing matrix
w = spalloc(3*n, 3*n, 27*n);

% Populate smoothing matrix
for i = 1:n % For each element, 
   for j = -2:0 % For each of the 3 components,
      w(3*i+j, 3*i+j) = 3; % Diagonal components = 3
      w(3*i+j, 3*share(i, share(i, :) ~= 0)+j) = -1; % Off diagonal components = -1
   end
end   
