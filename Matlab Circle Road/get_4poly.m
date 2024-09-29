function poly = get_4poly(Vec)

p = 3;
N = length(Vec);

combs=nchoosek(1:N,p);
m=size(combs,1);
S=sparse([1:m;1:m],combs.',1,m,N);
exp_comb = [S;speye(N);2*speye(N)];

poly = prod(Vec.^exp_comb, 2);
end