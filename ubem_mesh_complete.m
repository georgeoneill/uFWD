function surf = ubem_mesh_complete(meshes)

% Function to do some bookeeping on the previously generated surfaces,
% assumes that these are correct (ie normals point outwards, meshes have a
% solid angle of 4pi sr etc.)


type = {'brain','skull','scalp'};


for ii = 1:length(meshes)
    
    surf(ii).name           = type(ii);
    surf(ii).pos            = meshes(ii).pos;
    surf(ii).npos           = length(meshes(ii).pos);
    surf(ii).tri            = meshes(ii).tri;
    surf(ii).ntri           = length(meshes(ii).tri);

    
    pos     = meshes(ii).pos;
    tri    = meshes(ii).tri;
    p1      = pos(tri(:,1),:);
    p2      = pos(tri(:,2),:);
    p3      = pos(tri(:,3),:);
    
    surf(ii).tri_cent       = (p1+p2+p3)/3;
    tmp                     = cross((p2-p1),(p3-p1));
    surf(ii).tri_area       = 0.5*(vnorm(tmp,2));
    surf(ii).tri_nrms       = tmp./(2*surf(ii).tri_area*ones(1,3));
    surf(ii).tri_neighbours = triangle_neighbours(surf(ii).tri,surf(ii).npos);
    surf(ii).nrms           = normstri2vert(surf(ii).tri,surf(ii).tri_nrms,surf(ii).tri_area,surf(ii).npos);
    
    
end


end

function neighbour_tri = triangle_neighbours(tri,npos)
sz = size(tri);
neighbour_tri = cell(1,npos);
for ii = 1:npos
    [idx, ~] = ind2sub(sz,find(tri==ii));
    neighbour_tri{ii} = idx;
end

end

function vnorms = normstri2vert(tri,trinorms,triarea,npos)

verts = reshape(tri,numel(tri),1);
s = size(tri);
vnorms = zeros(3,npos);
for ii = 1:npos
    [idx, ~] = ind2sub(s,find(verts==ii));
    vnorms(:,ii) = sum(trinorms(idx,:).*(triarea(idx)*ones(1,3)));
end

vnorms = (vnorms./(ones(3,1)*vnorm(vnorms,1)))';

end

function size = vnorm(vec,dim)
size = sqrt(sum(vec.^2,dim));
end