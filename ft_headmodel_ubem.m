function headmodel = ft_headmodel_ubem(meshes, varargin)

% Supply conductivity values depending on number of meshes.
switch length(meshes)
    case 3
        opts.cond = ft_getopt(varargin,'conductivity',[0.3, 0.006, 0.3]);
    case 1
        opts.cond = ft_getopt(varargin,'conductivity',0.3);
end

% Just checks the meshes use tri/pos.
meshes = ubem_mesh_rename_assets(meshes);

% Check the meshes are in m units for now, life will be *so* much simpler
% later!
meshes = ft_determine_units(meshes);
meshes = ft_convert_units(meshes,'m');

% Check mesh(es) to see that they are good to go.
% Tests will include:
%               1) does the solid angle of each mesh come to +4pi sr?
%               2) is mesh 1 inside 2, which are both inside 3?
%               3) more tests can come.
opts.check_meshes = ft_getopt(varargin,'check_meshes',1);
if opts.check_meshes
    for ii = 1:length(meshes)
        try
            [~, meshes(ii)] = ubem_mesh_check(meshes(ii),1);
        catch
            % Will try to repair meshes in case of failure
            fprintf('Initial sanity check of mesh %d failed, repairing...',ii)
            if isempty(which('meshcheckrepair'))
                ft_hastoolbox('iso2mesh',1)
            end
            [meshes(ii).vertices, meshes(ii).faces] = meshcheckrepair(meshes(ii).vertices, meshes(ii).faces);
            disp('checking again')
            ubem_util_check_meshes(meshes(ii),1);
        end
    end
    
    if length(meshes) > 1
        % check meshes are in the correct order, if not, flip!
        try
            ubem_mesh_check(meshes,2);
        catch
            disp('meshes may be nested in the wrong order, reversing...')
            tmp = meshes;
            meshes(1) = tmp(3);
            meshes(2) = tmp(2);
            meshes(3) = tmp(1);
            ubem_mesh_check(meshes,2);
        end
    end
    
end

% Add addtional geometric information about meshes, and then rename
% everyhing, again. (This is fairly lightweight so can be done every time
% we need it).
surf = ubem_mesh_complete(meshes);

% Solve BEM for M/EEG - Adding sensors or electrodes come later.
bem = ubem_solve_bem(surf,opts.cond);

% Pack up
headmodel       = [];
headmodel.bnd   = meshes; % bnd order modified if they were the wrong way around
headmodel.type  = 'ubem';
headmodel.ubem  = bem;
headmodel       = ft_determine_units(headmodel);
