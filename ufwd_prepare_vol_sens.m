function [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, varargin)

% FT_PREPARE_VOL_SENS does some bookkeeping to ensure that the volume
% conductor model and the sensor array are ready for subsequent forward
% leadfield computations. It takes care of some pre-computations that can
% be done efficiently prior to the leadfield calculations.
%
% Use as
%   [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, ...)
% with input arguments
%   headmodel = structure with volume conductor definition
%   sens      = structure with gradiometer or electrode definition
%
% The headmodel structure represents a volume conductor model of the head,
% its contents depend on the type of model. It is described in more detail
% in FT_DATATYPE_HEADMODEL. The sens structure represents a electrode or
% gradiometer array. It is described in more detail in FT_DATATYPE_SENS.
%
% Additional options should be specified in key-value pairs and can be
%   'channel'  = cell-array with strings (default = 'all')
%
% The detailed behavior of this function depends on whether the input
% consists of EEG or MEG and furthermoree depends on the type of volume
% conductor model:
% - in case of EEG single and concentric sphere models, the electrodes are
%   projected onto the skin surface.
% - in case of EEG boundary element models, the electrodes are projected on
%   the surface and a blilinear interpoaltion matrix from vertices to
%   electrodes is computed.
% - in case of MEG and a localspheres model, a local sphere is determined
%   for each coil in the gradiometer definition.
%  - in case of MEG with a singleshell Nolte model, the volume conduction
%    model is initialized
% In any case channel selection and reordering will be done. The channel
% order returned by this function corresponds to the order in the 'channel'
% option, or if not specified, to the order in the input sensor array.
%
% See also FT_COMPUTE_LEADFIELD, FT_READ_HEADMODEL, FT_READ_SENS

% Copyright (C) 2004-2015, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if iscell(headmodel) && iscell(sens)
    % this represents combined EEG, ECoG and/or MEG
    for i=1:numel(headmodel)
        [headmodel{i}, sens{i}] = ft_prepare_vol_sens(headmodel{i}, sens{i}, varargin{:});
    end
    return
end

% get the optional input arguments
% fileformat = ft_getopt(varargin, 'fileformat');
channel = ft_getopt(varargin, 'channel', sens.label);   % cell-array with channel labels, default is all

% ensure that the sensor description is up-to-date (Aug 2011)
sens = ft_datatype_sens(sens);

% this is to support volumes saved in mat-files, particularly interpolated
if ischar(headmodel)
    vpath     = fileparts(headmodel);   % remember the path to the file
    headmodel = ft_read_headmodel(headmodel); % replace the filename with the content of the file
end

% ensure that the volume conduction description is up-to-date (Jul 2012)
headmodel = ft_datatype_headmodel(headmodel);

% determine whether the input contains EEG or MEG sensors
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

% determine the skin compartment
if ~isfield(headmodel, 'skin_surface')
    if isfield(headmodel, 'bnd')
        headmodel.skin_surface   = find_outermost_boundary(headmodel.bnd);
    elseif isfield(headmodel, 'r') && length(headmodel.r)<=4
        [dum, headmodel.skin_surface] = max(headmodel.r);
    end
end

% determine the inner_skull_surface compartment
if ~isfield(headmodel, 'inner_skull_surface')
    if isfield(headmodel, 'bnd')
        headmodel.inner_skull_surface  = find_innermost_boundary(headmodel.bnd);
    elseif isfield(headmodel, 'r') && length(headmodel.r)<=4
        [dum, headmodel.inner_skull_surface] = min(headmodel.r);
    end
end

% FT_HEADMODELTYPE to an empty struct won't work further down
if isempty(headmodel)
    headmodel = [];
end

% this makes them easier to recognise
sens.type       = ft_senstype(sens);
headmodel.type  = ft_headmodeltype(headmodel);

if isfield(headmodel, 'unit') && isfield(sens, 'unit') && ~strcmp(headmodel.unit, sens.unit)
    ft_error('inconsistency in the units of the volume conductor and the sensor array');
end

if ismeg && iseeg
    % this is something that could be implemented relatively easily
    ft_error('simultaneous EEG and MEG not yet supported');
    
elseif ~ismeg && ~iseeg
    ft_error('the input does not look like EEG, nor like MEG');
    
elseif ismeg
    
    % always ensure that there is a linear transfer matrix for combining the coils into gradiometers
    if ~isfield(sens, 'tra');
        Nchans = length(sens.label);
        Ncoils = size(sens.coilpos,1);
        if Nchans~=Ncoils
            ft_error('inconsistent number of channels and coils');
        end
        sens.tra = eye(Nchans, Ncoils);
    end
    
    if ~ft_headmodeltype(headmodel, 'localspheres')
        % select the desired channels from the gradiometer array
        [selchan, selsens] = match_str(channel, sens.label);
        % only keep the desired channels, order them according to the users specification
        try, sens.chantype = sens.chantype(selsens,:); end
        try, sens.chanunit = sens.chanunit(selsens,:); end
        try, sens.chanpos  = sens.chanpos (selsens,:); end
        try, sens.chanori  = sens.chanori (selsens,:); end
        sens.label    = sens.label(selsens);
        sens.tra      = sens.tra(selsens,:);
    else
        % for the localspheres model it is done further down
    end
    
    % remove the coils that do not contribute to any channel output
    selcoil      = any(sens.tra~=0,1);
    sens.coilpos = sens.coilpos(selcoil,:);
    sens.coilori = sens.coilori(selcoil,:);
    sens.tra     = sens.tra(:,selcoil);
    
    switch ufwd_headmodeltype(headmodel)
        
        case 'ubem'
            % Generate the transfer matrix based on the sensor information
            sens    = ft_convert_units(sens,headmodel.unit);
            coils.r = sens.coilpos;
            coils.o = sens.coilori;
            % need to repopulate the additional geometry information for
            % the meshes then compute transfer matrix
            surf    = ubem_mesh_complete(headmodel.bnd);
            xMEG    = ubem_solve_coils(surf,headmodel.ubem,coils);
            % housekeeping - we can save memory by only keeping the
            % essentials, which is just the transfer matrix
            headmodel               = rmfield(headmodel,'ubem');
            headmodel.ubem.xMEG     = xMEG;
            headmodel.ubem.surfs    = surf;
            
    end
    
elseif iseeg
    
    % the electrodes are used, the channel positions are not relevant any more
    % channel positinos need to be recomputed after projecting the electrodes on the skin
    if isfield(sens, 'chanpos'); sens = rmfield(sens, 'chanpos'); end
    
    % select the desired channels from the electrode array
    % order them according to the users specification
    [selchan, selsens] = match_str(channel, sens.label);
    Nchans = length(sens.label);
    
    sens.label     = sens.label(selsens);
    try, sens.chantype  = sens.chantype(selsens); end
    try, sens.chanunit  = sens.chanunit(selsens); end
    
    if isfield(sens, 'tra')
        % first only modify the linear combination of electrodes into channels
        sens.tra     = sens.tra(selsens,:);
        % subsequently remove the electrodes that do not contribute to any channel output
        selelec      = any(sens.tra~=0,1);
        sens.elecpos = sens.elecpos(selelec,:);
        sens.tra     = sens.tra(:,selelec);
    else
        % the electrodes and channels are identical
        sens.elecpos = sens.elecpos(selsens,:);
    end
    
    switch ufwd_headmodeltype(headmodel)
        
        case 'ubem'
            
            els.r = sens.elecpos;      
            % need to repopulate the additional geometry information for
            % the meshes then compute transfer matrix
            surf    = ubem_mesh_complete(headmodel.bnd);
            [xEEG, prj]    = ubem_solve_electrodes(surf,headmodel.ubem,els);
            % housekeeping - we can save memory by only keeping the
            % essentials, which is just the transfer matrix
            headmodel               = rmfield(headmodel,'ubem');
            headmodel.ubem.xEEG     = xEEG;
            headmodel.ubem.surfs    = surf;
            sens.elecpos            = prj;     
              
        case {'bem', 'dipoli', 'asa', 'bemcp'}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % do postprocessing of volume and electrodes in case of BEM model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % project the electrodes on the skin and determine the bilinear interpolation matrix
            if ~isfield(headmodel, 'tra') && (isfield(headmodel, 'mat') && ~isempty(headmodel.mat))
                % determine boundary corresponding with skin and inner_skull_surface
                if ~isfield(headmodel, 'skin_surface')
                    headmodel.skin_surface = find_outermost_boundary(headmodel.bnd);
                    fprintf('determining skin compartment (%d)\n', headmodel.skin_surface);
                end
                if ~isfield(headmodel, 'source')
                    headmodel.source = find_innermost_boundary(headmodel.bnd);
                    fprintf('determining source compartment (%d)\n', headmodel.source);
                end
                if size(headmodel.mat,1)~=size(headmodel.mat,2) && size(headmodel.mat,1)==length(sens.elecpos)
                    fprintf('electrode transfer and system matrix were already combined\n');
                else
                    fprintf('projecting electrodes on skin surface\n');
                    % compute linear interpolation from triangle vertices towards electrodes
                    [el, prj] = project_elec(sens.elecpos, headmodel.bnd(headmodel.skin_surface).pos, headmodel.bnd(headmodel.skin_surface).tri);
                    tra       = transfer_elec(headmodel.bnd(headmodel.skin_surface).pos, headmodel.bnd(headmodel.skin_surface).tri, el);
                    
                    % replace the original electrode positions by the projected positions
                    sens.elecpos = prj;
                    
                    if size(headmodel.mat,1)==size(headmodel.bnd(headmodel.skin_surface).pos,1)
                        % construct the transfer from only the skin vertices towards electrodes
                        interp = tra;
                    else
                        % construct the transfer from all vertices (also inner_skull_surface/outer_skull_surface) towards electrodes
                        interp = [];
                        for i=1:length(headmodel.bnd)
                            if i==headmodel.skin_surface
                                interp = [interp, tra];
                            else
                                interp = [interp, zeros(size(el,1), size(headmodel.bnd(i).pos,1))];
                            end
                        end
                    end
                    
                    % incorporate the linear interpolation matrix and the system matrix into one matrix
                    % this speeds up the subsequent repeated leadfield computations
                    fprintf('combining electrode transfer and system matrix\n');
                    
                    % convert to sparse matrix to speed up the subsequent multiplication
                    interp  = sparse(interp);
                    headmodel.mat = interp * headmodel.mat;
                    % ensure that the model potential will be average referenced
                    avg = mean(headmodel.mat, 1);
                    headmodel.mat = headmodel.mat - repmat(avg, size(headmodel.mat,1), 1);
                end
            end
            
            keyboard
            
        case  'openmeeg'
            % don't do anything, h2em or h2mm generated later in ft_prepare_leadfield
            
        case 'fns'
            if isfield(headmodel,'bnd')
                [el, prj] = project_elec(sens.elecpos, headmodel.bnd.pos, headmodel.bnd.tri);
                sens.tra = transfer_elec(headmodel.bnd.pos, headmodel.bnd.tri, el);
                % replace the original electrode positions by the projected positions
                sens.elecpos = prj;
            end
            
        case 'simbio'
            % check that the external toolbox is present
            ft_hastoolbox('simbio', 1);
            
            % extract the outer surface
            bnd = mesh2edge(headmodel);
            for j=1:length(sens.label)
                d = bsxfun(@minus, bnd.pos, sens.elecpos(j,:));
                [d, i] = min(sum(d.^2, 2));
                % replace the position of each electrode by the closest vertex
                sens.elecpos(j,:) = bnd.pos(i,:);
            end
            
            if (isfield(headmodel,'transfer') && isfield(headmodel,'elec'))
                if all(ismember(sens.label,headmodel.elec.label))
                    [sensmember, senslocation] = ismember(sens.label,headmodel.elec.label);
                    if (norm(sens.elecpos - headmodel.elec.elecpos(senslocation,:))<1e-8)
                        headmodel.transfer = headmodel.transfer(senslocation,:);
                        headmodel.elec = sens;
                    else
                        ft_error('Electrode positions do not fit to the given transfer matrix!');
                    end
                else
                    ft_error('Transfer matrix does not fit the given set of electrodes!');
                end
            else
                headmodel.transfer = sb_transfer(headmodel,sens);
                headmodel.elec = sens;
            end
            
            
        case 'interpolate'
            % this is to allow moving leadfield files
            if ~exist(headmodel.filename{1}, 'file')
                for i = 1:length(headmodel.filename)
                    [p, f, x] = fileparts(headmodel.filename{i});
                    headmodel.filename{i} = fullfile(vpath, [f x]);
                end
            end
            
            matchlab = isequal(sens.label, headmodel.sens.label);
            matchpos = isequal(sens.elecpos, headmodel.sens.elecpos);
            matchtra = (~isfield(sens, 'tra') && ~isfield(headmodel.sens, 'tra')) || isequal(sens.tra, headmodel.sens.tra);
            
            if matchlab && matchpos && matchtra
                % the input sensor array matches precisely with the forward model
                % no further interpolation is needed
            else
                % interpolate the channels in the forward model to the desired channels
                filename = tempname;
                headmodel  = ft_headmodel_interpolate(filename, sens, headmodel);
                % update the sensor array with the one from the volume conductor
                sens = headmodel.sens;
            end % if recomputing interpolation
            
            % for the leadfield computations the @nifti object is used to map the image data into memory
            ft_hastoolbox('spm8up', 1);
            for i=1:length(headmodel.sens.label)
                % map each of the leadfield files into memory
                headmodel.chan{i} = nifti(headmodel.filename{i});
            end
            
        otherwise
            ft_error('unsupported volume conductor model for EEG');
    end
    
    % FIXME this needs careful thought to ensure that the average referencing which is now done here and there, and that the linear interpolation in case of BEM are all dealt with consistently
    % % always ensure that there is a linear transfer matrix for
    % % rereferencing the EEG potential
    % if ~isfield(sens, 'tra');
    %   sens.tra = eye(length(sens.label));
    % end
    
    % update the channel positions as the electrodes were projected to the skin surface
    [pos, ori, lab] = channelposition(sens);
    [selsens, selpos] = match_str(sens.label, lab);
    sens.chanpos = nan(length(sens.label),3);
    sens.chanpos(selsens,:) = pos(selpos,:);
    
end % if iseeg or ismeg

if isfield(sens, 'tra')
    if issparse(sens.tra) && size(sens.tra, 1)==1
        % this multiplication would result in a sparse leadfield, which is not what we want
        % the effect can be demonstrated as sparse(1)*rand(1,10), see also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1169#c7
        sens.tra = full(sens.tra);
    elseif ~issparse(sens.tra) && size(sens.tra, 1)>1
        % the multiplication of the "sensor" leadfield (electrode or coil) with the tra matrix to get the "channel" leadfield
        % is faster for most cases if the pre-multiplying weighting matrix is made sparse
        sens.tra = sparse(sens.tra);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ppr = pointproj(P,plane)
% projects a point on a plane
% plane(1:3) is a point on the plane
% plane(4:6) is the ori of the plane
Ppr  = [];
ori  = plane(4:6);
line = [P ori];
% get indices of line and plane which are parallel
par = abs(dot(plane(4:6), line(:,4:6), 2))<1e-14;
% difference between origins of plane and line
dp = plane(1:3) - line(:, 1:3);
% Divide only for non parallel vectors (DL)
t = dot(ori(~par,:), dp(~par,:), 2)./dot(ori(~par,:), line(~par,4:6), 2);
% compute coord of intersection point
Ppr(~par, :) = line(~par,1:3) + repmat(t,1,3).*line(~par,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as a replacement for the dist function in the Neural
% Networks toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d] = dist(x)
n = size(x,2);
d = zeros(n,n);
for i=1:n
    for j=(i+1):n
        d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
        d(j,i) = d(i,j);
    end
end
