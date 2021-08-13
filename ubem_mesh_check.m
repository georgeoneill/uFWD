function [pass, meshes]= ubem_mesh_check(meshes,test)

if nargin==1
    test=1;
end

switch test
    
    case 0 % Silent running version of test 1, are triangles correctly oriented?
        
        for ii = 1:length(meshes)
            sr = get_solids_fast(meshes(ii).pos,meshes(ii).tri);
            if abs(sr - (2*pi)) < 1e-6
                pass(ii) = 1;
            elseif abs(sr + (2*pi)) < 1e-6
                pass(ii) = -1;
            else
                pass(ii) = 0;
            end
        end
        
    case 1
        
        % Test 1: are the surfaces complete?
        fprintf('Checking boundary mesh integrity: ')
        tmp_pass = ubem_mesh_check(meshes,0);
        if tmp_pass == 1
            fprintf('OK!\n')
            pass = 1;
        elseif tmp_pass == -1
            fprintf('inside out: correcting: ')
            ftmp = meshes.tri;
            meshes.tri = ftmp(:,[1 3 2]);
            tmp_pass_2 = ubem_mesh_check(meshes,0);
            if tmp_pass_2~=1
                fprintf('FAIL\n')
                error('Boundaries are not completely closed surfaces!');
            else
                fprintf('OK!\n')
                pass = 1;
            end
        else
            fprintf('FAIL\n')
            error('Boundaries are not completely closed surfaces!');
        end
        
    case 2
        
        % Test 2: Is surface 1 entirely within 2, which itself is in 3?
        % Obvs, skip the test is there is a single layer BEM.
        if length(meshes) > 1
            fprintf('Checking boundaries are entirely nested within each other: ')
            clear tmp_pass
            for ii = 1:length(meshes)-1
                fro = meshes(ii+1);
                to = meshes(ii);
                clear sr
                sr = get_solids_fast(fro.pos,fro.tri,to.pos);
                in = (abs(sr - (2*pi)) < 1e-4);
                tmp_pass(ii) = sum(in) == numel(to.pos)./3;
            end
            pass = sum(tmp_pass) == length(meshes)-1;
            if pass
                fprintf('OK!\n')
            else
                fprintf('FAIL\n')
                error('One or more are boundaries are not nested properly');
            end
        else
            pass = 1;
        end
        
end

end

function [solids] = get_solids_fast(pos,tri,varargin)

if nargin==3
    fros = varargin{1};
    npoints = length(fros);
else
    fros = mean(pos,1);
    npoints = 1;
end

% need to order the verices in order of all the triangles.

vtris = zeros(length(tri),3,3);
for ii = 1:length(tri)
    vtris(ii,:,:) = pos(tri(ii,:),:);
end

slices = arange(1,100,npoints);
if npoints == 1
    slices = [1 1];
end

solids = zeros(1,npoints);
for ii = 1:length(slices)-1
    tmp1 = reshape(fros(slices(ii):slices(ii+1),:),1,1,numel(slices(ii):slices(ii+1)),[]);
    tmp2 = reshape(permute(vtris,[2 1 3]),size(vtris,2),size(vtris,1),1,size(vtris,3));
    vs = bsxfun(@minus,tmp1,tmp2);
    v1 = squeeze(vs(1,:,:,:));
    v2 = squeeze(vs(2,:,:,:));
    v3 = squeeze(vs(3,:,:,:));
    triples = fast_cross_nd_sum(v1,v2,v3);
    ls = vnorm(vs,4);
    ss = squeeze(prod(ls,1));
    switch(ndims(v1))
        case 2 % just incase there is one point to compare to.
            s1 = ls(3,:,:)'.*sum(v1.*v2,ndims(v1));
            s2 = ls(2,:,:)'.*sum(v1.*v3,ndims(v1));
            s3 = ls(1,:,:)'.*sum(v2.*v3,ndims(v1));
            ss = ss'+s1+s2+s3;
            solids(slices(ii):slices(ii+1)) = -sum(bsxfun(@atan2,triples,ss));
        case 3
            s1 = squeeze(ls(3,:,:)).*sum(v1.*v2,ndims(v1));
            s2 = squeeze(ls(2,:,:)).*sum(v1.*v3,ndims(v1));
            s3 = squeeze(ls(1,:,:)).*sum(v2.*v3,ndims(v1));
            ss = ss+s1+s2+s3;
            solids(slices(ii):slices(ii+1)) = -sum(bsxfun(@atan2,triples,squeeze(ss)),1);
    end
end

end

function X = fast_cross(a,b)
i = a(:,2).*b(:,3)-a(:,3).*b(:,2);
j = a(:,3).*b(:,1)-a(:,1).*b(:,3);
k = a(:,1).*b(:,2)-a(:,2).*b(:,1);
X = [i j k];
end

function triples = fast_cross_nd_sum(a,b,c)
otherdims = repmat({':'},1,ndims(a)-1);
triples = (a(otherdims{:},2).*b(otherdims{:},3)-a(otherdims{:},3).*b(otherdims{:},2)).*c(otherdims{:},1)...
    + (a(otherdims{:},3).*b(otherdims{:},1)-a(otherdims{:},1).*b(otherdims{:},3)).*c(otherdims{:},2)...
    + (a(otherdims{:},1).*b(otherdims{:},2)-a(otherdims{:},2).*b(otherdims{:},1)).*c(otherdims{:},3);
end

function vec = arange(lo,step,hi)
vec = lo:step:hi;
if vec(end) ~= hi
    vec = [vec hi];
end
end

function size = vnorm(vec,dim)
size = sqrt(sum(vec.^2,dim));
end