function meshes2 = ubem_mesh_rename_assets(meshes)

for ii = 1:length(meshes)
    
    if isfield(meshes(ii),'pos'); meshes2(ii).pos = meshes(ii).pos; end
    if isfield(meshes(ii),'pnt'); meshes2(ii).pos = meshes(ii).pnt; end
    if isfield(meshes(ii),'vertices'); meshes2(ii).pos = meshes(ii).vertices; end
    if isfield(meshes(ii),'tri'); meshes2(ii).tri = meshes(ii).tri; end
    if isfield(meshes(ii),'faces'); meshes2(ii).tri = meshes(ii).faces; end
    if isfield(meshes(ii),'name'); meshes2(ii).name = meshes(ii).name; end
    if isfield(meshes(ii),'unit'); meshes2(ii).unit = meshes(ii).unit; end
    
end