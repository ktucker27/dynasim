function sa = get_sa(s, n, a)

%sa = kron(kron(speye(2^(n-a), 2^(n-a)),s), speye(2^(a-1), 2^(a-1)));
sa = kron(kron(speye(2^(a-1), 2^(a-1)),s), speye(2^(n-a), 2^(n-a)));