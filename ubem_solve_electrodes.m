function [xEEG, prj] = ubem_solve_electrodes(surf,bem,els,varargin)

% electrode transfer matrix generator, using the projection algorthms from
% FieldTrip!

fprintf('Generating EEG transfer matrix: ')

if nargin == 4
    skin_id = varargin{1};
else
    skin_id = 3; 
end

[el, prj]       = project_elec(els.r, surf(skin_id).pos, surf(skin_id).tri);
tra             = transfer_elec(surf(skin_id).pos, surf(skin_id).tri, el);

% construct the interpolator
interp = [];
for ii=1:length(surf)
    if ii==skin_id
        interp = [interp, tra];
    else
        interp = [interp, zeros(size(el,1), surf(ii).npos)];
    end
end

% collapse down the transfer matrix
xEEG = interp*bem.solution;

% scale
mults = [];
for ii = 1:length(surf)
    mults = cat(1,mults,bem.source_mult(ii)*ones(surf(ii).npos,1));
end
mults = mults/(4*pi);
xEEG = bsxfun(@times,xEEG,mults');

fprintf('COMPLETE\n')