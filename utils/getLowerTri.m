function [list idx] = getLowerTri(rdm)

[s1 s2] = size(rdm);
assert(s1==s2)

idx = find(tril(ones(s1), -1));

list = rdm(idx);
