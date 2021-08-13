function bem = ubem_solve_bem(surf,conductivity)

bem = struct;
bem.sigma = conductivity; %Condutvitity from outside to inside.
bem = add_gamma_multipliers(bem);
coeffs = bem_lin_pot_coeff(surf);
disp('Inverting Solution')
bem.solution = bem_multi_solution(coeffs, bem.gamma, [surf.npos]);
if length(surf) == 3 &&  bem.sigma(2)/bem.sigma(1) < 0.1
    % Take the innermost compartment, recalculate coefficients and modify
    % the original solution accordingly.
    disp('IP Approach required')
    coeffs = bem_lin_pot_coeff(surf(1));
    disp('Inverting inner most shell, again')
    IP = bem_multi_solution(coeffs, [], surf(1).npos);
    disp('Modifying solution')
    bem.solution = bem_ip_modify_solution(bem.solution,IP,bem.sigma(2)/bem.sigma(1),[surf.npos]);
end
disp('BEM Solved')

end

function coeff = bem_lin_pot_coeff(surf)

% Calculate the coefficients for linear collocation approach.
nsurfs = length(surf);
npos = [surf.npos];
coeff = zeros(sum([surf.npos]));
offsets = [0 cumsum(npos)];
fprintf('calculating potential coefficients:\n')
for ii = 1:nsurfs
    
    pos_ord = 1:npos(ii);
    
    for jj = 1:nsurfs
        
        tmp = zeros(surf(ii).npos,surf(jj).npos);
        fprintf('\t%s -> %s\n',cell2mat(surf(ii).name),cell2mat(surf(jj).name));
        tri_nrms    = surf(jj).tri_nrms;
        tri_area    = surf(jj).tri_area;
        
        for kk = 1:surf(jj).ntri
            tri         = surf(jj).tri(kk,:);
            tri_pos      = surf(jj).pos(tri,:);
            if ii == jj
                skip_idx = pos_ord == tri(1) | pos_ord == tri(2) | pos_ord == tri(3);
            else 
                skip_idx = [];
            end
            coeffs = lin_pot_coeff(surf(ii).pos, tri_pos, tri_nrms(kk,:), tri_area(kk));
            coeffs(skip_idx,:) = 0;
            tmp(:,tri) = tmp(:,tri) - coeffs;
        end
        
        if ii == jj
            tmp = correct_auto_elements(surf(ii),tmp);
        end
        
        % check for NaNs or Infs.
        chk = isnan(tmp) | isinf(tmp);
        if sum(chk(:)) > 0
            error('NaNs or Infs detected in potential coefficients, check meshes!')
        end
        
        coeff((offsets(ii)+1):(offsets(ii+1)),(offsets(jj)+1):(offsets(jj+1))) = tmp;
        
    end

end
end

function solids_inv = bem_multi_solution(solids, gamma, nps)

nsurf = length(nps);
pi2 = 1/(2*pi);
n_tot = sum(nps);
defl = 1/n_tot;
offsets = [0 cumsum(nps)];

tmp = zeros(size(solids));

for ii = 1:nsurf
    for jj = 1:nsurf
        
        if ~isempty(gamma)
            mult = pi2*gamma(ii,jj);
        else
            mult = pi2;
        end
        tmp((offsets(ii)+1):(offsets(ii+1)),(offsets(jj)+1):(offsets(jj+1))) =...
            defl - mult*solids((offsets(ii)+1):(offsets(ii+1)),(offsets(jj)+1):(offsets(jj+1)));
    end
end
   
tmp = tmp+eye(size(tmp));
solids_inv = inv(tmp);


end

function fixed = bem_ip_modify_solution(broken,ip,ip_mult,nps)


n_last = nps(end);
offsets = [0 cumsum(nps)];
mult = (1+ip_mult)/ip_mult;
disp('Combining')
for ii = 1:length(nps)
    % The BEM solution can be thought of as 3x3 block matrix, the left
    % most column of the matrix needs modifying block by block
    
    tmp = broken((offsets(ii)+1):(offsets(ii+1)),(offsets(1)+1):(offsets(2)));
    tmp = tmp - 2*tmp*ip;
    broken((offsets(ii)+1):(offsets(ii+1)),(offsets(1)+1):(offsets(2))) = tmp;
        
end

% special treatment for the top left block
tmp = broken((offsets(1)+1):(offsets(2)),(offsets(1)+1):(offsets(2)));
tmp = tmp + mult*ip;
broken((offsets(1)+1):(offsets(2)),(offsets(1)+1):(offsets(2))) = tmp;

fixed = ip_mult*broken;

end

function omega = lin_pot_coeff(pos, tri_pos, tri_nrms, tri_area)

omega = zeros(size(pos));

v1 = ones(length(pos),1)*tri_pos(1,:) - pos;
v2 = ones(length(pos),1)*tri_pos(2,:) - pos;
v3 = ones(length(pos),1)*tri_pos(3,:) - pos;
triples = fast_cross_nd_sum(v1,v2,v3);
l1 = vnorm(v1,2);
l2 = vnorm(v2,2);
l3 = vnorm(v3,2);
ss = l1.*l2.*l3;
s1 = l3.*sum(v1.*v2,2);
s2 = l2.*sum(v1.*v3,2);
s3 = l1.*sum(v2.*v3,2);
ss = ss+s1+s2+s3;
solids = bsxfun(@atan2,triples,ss);
bad_mask = (abs(solids) < (pi/1e6));
l1(bad_mask) = 1;
l2(bad_mask) = 1;
l3(bad_mask) = 1;

beta = zeros(length(pos),3);
beta(:,1) = calc_beta(v1,l1,v2,l2);
beta(:,2) = calc_beta(v2,l2,v3,l3);
beta(:,3) = calc_beta(v3,l3,v1,l1);


vec_omega = ((beta(:,3) - beta(:,1))*ones(1,3)).*v1;
vec_omega = vec_omega + ((beta(:,1) - beta(:,2))*ones(1,3)).*v2;
vec_omega = vec_omega + ((beta(:,2) - beta(:,3))*ones(1,3)).*v3;

area2 = 2*tri_area;
n2 = 1/(area2.*area2);
yys = cat(3,v1,v2,v3);
idx = [3 1 2; 2 3 1];

for ii = 1:3
    diff = yys(:,:,idx(1,ii)) - yys(:,:,idx(2,ii));
    zdots = fast_cross_nd_sum(yys(:,:,idx(2,ii)),yys(:,:,idx(1,ii)),tri_nrms);
    omega(:,ii) = -n2 * ( 2 * area2 * zdots .* solids - triples .* sum(diff.*vec_omega,2));
end
omega(bad_mask,:) = 0;

end

function mat = correct_auto_elements(surf,mat)
pi2 = 2*pi;
misses = pi2 - sum(mat,2);
for ii = 1:length(misses)
    miss = misses(ii);
    n_memb = length(surf.tri_neighbours{ii});
    mat(ii,ii) = miss/2;
    miss = miss/(2*n_memb);
    neighbours = unique(surf.tri(surf.tri_neighbours{ii},:));
    neighbours(neighbours==ii) = [];
    for jj = 1:length(neighbours)   
        mat(ii,neighbours(jj)) = mat(ii,neighbours(jj))+miss;
    end
end
end

function triples = fast_cross_nd_sum(a,b,c)
    triples = (a(:,2).*b(:,3)-a(:,3).*b(:,2)).*c(:,1)...
        + (a(:,3).*b(:,1)-a(:,1).*b(:,3)).*c(:,2)...
        + (a(:,1).*b(:,2)-a(:,2).*b(:,1)).*c(:,3);
end

function beta = calc_beta(p1,p1_norm,p2,p2_norm)
p21    = p2(1,:) - p1(1,:);
sz      = vnorm(p21,2);
p21    = p21./sz;
tmp     = sum(bsxfun(@times,[p1;p2],p21),2);
num     = p1_norm + tmp(1:length(p1));
den     = p2_norm + tmp(length(p1)+1:end);
tmp1    = num./den;
tmp2    = log(tmp1);
beta    = tmp2./sz;
end

function bem = add_gamma_multipliers(bem)

if length(bem.sigma) == 3
    sigma = [bem.sigma 0];
    bem.source_mult = 2./(sigma(1:3)+sigma(2:end));
    bem.field_mult = sigma(1:3)-sigma(2:end);
    num = (sigma(1:3)-sigma(2:end));
    den = (sigma(1:3)+sigma(2:end));
    bem.gamma = bsxfun(@rdivide,num,den');
elseif length(bem.sigma) == 1
    bem.source_mult = 2./bem.sigma;
    bem.field_mult = bem.sigma;
    bem.gamma = 1;
else
    error('meh')
end

end

function size = vnorm(vec,dim)
size = sqrt(sum(vec.^2,dim));
end
