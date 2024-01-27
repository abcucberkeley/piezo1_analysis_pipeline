function GaussInt = GU_calcGaussianIntegral(A, s, varargin)
% calculated gaussian integral
% input: Amplitude (single or an array); vector of sigmas;
% n-d requires n sigma input
% s must be in rows for more than 1 value
% Gokul Upadhyayula, March, 2017


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('A', @isvector);
ip.addRequired('s', @ismatrix);
ip.parse(A, s, varargin{:});

nd = numel(s)/numel(A);
GaussInt = A.*((2*pi)^(nd/2)* prod(s));

% %  test code below
% 
% opts.Alpha = 0.05;
% opts.CellMask = [];
% opts.RemoveRedundant = true;
% opts.WindowSize = [];
% opts.Mode = 'xyzAsrc';
% 
% sig = 0.4;
% for i = 1:10
% a = zeros(51,51,51);
% a(25,25,25) = i;
% ag = filterGauss3D(a,sig);
% % imtool3D(ag)
% s(i) = max(ag(:));
% 
% [pstruct, mask] = pointSourceDetection3D(ag, 1, 'Mode', opts.Mode, 'Alpha', opts.Alpha,...
% 'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant,...
% 'RefineMaskLoG', false, 'WindowSize', opts.WindowSize);
% A(i) = pstruct.A;
% end
% GaussInt = GU_calcGaussianIntegral(A, repmat(sig,3,1));
