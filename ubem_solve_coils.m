function xMEG = ubem_solve_coils(surf,bem,coils)

fprintf('Generating MEG transfer matrix: ')

offsets = [0 cumsum([surf.npos])];
xMEG = zeros(length(coils.r),sum([surf.npos]));
for ii = 1:length(surf)
    coeff = lin_field_coeff(surf(ii),coils);
    coeff = coeff*bem.field_mult(ii);
    tmp = bem.solution((offsets(ii)+1):(offsets(ii+1)),:);
    xMEG = xMEG + coeff*tmp;
end

mults = [];
for ii = 1:length(surf)
    mults = cat(1,mults,bem.source_mult(ii)*ones(surf(ii).npos,1));
end
mults = mults/(4*pi);
xMEG = bsxfun(@times,xMEG,mults');

fprintf('COMPLETE\n')

end

function coeff = lin_field_coeff(surf,coils)

o = coils.o;
r = coils.r;

bem_pos = surf.pos;
tris = surf.tri;
tn = surf.tri_nrms;
ta = surf.tri_area;
coeff = zeros(length(r),length(bem_pos));

for ii = 1:surf.ntri
    tri = tris(ii,:);
    tri_p = bem_pos(tri,:);
    tnn = tn(ii,:);
    zz = zeros(length(r),3);
    for jj = 1:3
        tpp = tri_p(jj,:);
        diff = bsxfun(@minus,r,tpp);
        d2 = sum(diff.*diff,2);
        c = fast_cross(diff,tnn);
        x = ta(ii)*sum(c.*o,2)./(3.*d2.*sqrt(d2));
        zz(:,jj) = x;
    end
coeff(:,tri) = coeff(:,tri)+zz;

end
end

function c = fast_cross(a,b)

assert(size(a,2)==3)
assert(size(b,2)==3)

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) ...
    a(:,3).*b(:,1)-a(:,1).*b(:,3)...
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end
