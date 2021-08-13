function [bem, meshes] = ubem_headmodel(meshes, varargin)

% Supply conductivity values depending on number of meshes.
switch length(meshes)
    case 3
        opts.cond = ft_getopt(varargin,'conductivity',[0.3, 0.006, 0.3]);
    case 1
        opts.cond = ft_getopt(varargin,'conductivity',0.3);
end

% FieldTrip likes to use pos/pnt/tri for mashes rather than faces/vertices
% but all the code I have written expects the latter. This function renames
% all the relevent items in each mesh structre to the correct name
meshes = ubem_mesh_rename_assets(meshes);

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
            fprintf('Initial sanity check of %s mesh failed, repairing...',meshes(ii).name)
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

% With all sanity checks passed, lets fill in all the surface information
% based on just the faces and vertices alone. Then solve the geometric part
% of the BEM common to both EEG and MEG.
surf = ubem_mesh_complete(meshes);
bem = ubem_solve_bem(surf,opts.cond);
