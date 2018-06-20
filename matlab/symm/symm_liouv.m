function [L, Lfab, Sp, Sm, Spr, Smr] = symm_liouv(n, w, o, faa, fab, gaa, gab, gamma, gel)

d = (n+3)*(n+2)*(n+1)/6;

Lgaa = zeros(d,d);
Lgab = zeros(d,d);
Lfaa = zeros(d,d);
Lfab = zeros(d,d);
Lw = zeros(d,d);
Lel = zeros(d,d);

Sp = zeros(d,d);
Sm = zeros(d,d);
Spr = zeros(d,d);
Smr = zeros(d,d);

col_idx = [0;0;0];
for j=1:d
    nz = col_idx(1);
    np = col_idx(2);
    nm = col_idx(3);
    n1 = n - nz - np - nm;
    
    row_idx = [0;0;0];
    for i=1:d
        if idx_equals(row_idx, col_idx)
            Lw(i,j) = col_idx(1) + 0.5*(col_idx(2) + col_idx(3));
            %Lf(i,j) = 2*nz*nm + 2*nz + 2*nz*np + np + nm;
            Lfaa(i,j) = col_idx(1) + 0.5*(col_idx(2) + col_idx(3));
            Lfab(i,j) = 0.5*(nm - np)*n1 + 0.5*nz*(nm + 1 + 3*(np + 1)) ...
                + 0.5*(n1 + 1)*(np - nm) + 0.5*(nz+1)*(np + 3*nm);
            Lgaa(i,j) = col_idx(2) - col_idx(3);
            Lel(i,j) = col_idx(2) + col_idx(3);
        end
        
        if idx_equals(row_idx, col_idx + [1;0;0])
            Lgab(i,j) = -(col_idx(2) - col_idx(3))*(col_idx(1) + 1);
            Lw(i,j) = -(col_idx(1) + 1);
            %Lf(i,j) = nz*nm + 2*nz + nz*np + nm + 2 + np;
            Lfaa(i,j) = col_idx(1) + 1;
            Lfab(i,j) = 0.5*(nz+1)*(nm+1+3*(np+1)) - 0.5*(nz+1)*(np - nm);
        end
        
        if idx_equals(row_idx, col_idx + [-1;0;0])
            Lgab(i,j) = -(col_idx(2) - col_idx(3))*(n - sum(col_idx) + 1);
            Lfab(i,j) = -n1*nm - n1*np - nm - np;
            %Lfab(i,j) = 0.5*(n1+1)*(nm - np) + 0.5*(n1+1)*(nm-np);
        end
        
        if idx_equals(row_idx, col_idx + [-1;1;1])
            Lfab(i,j) = -4*(np+1)*(nm+1);
        end
        
        if idx_equals(row_idx, col_idx + [-2;1;1])
            Lfab(i,j) = -4*(np+1)*(nm+1);
        end
        
        if idx_equals(row_idx, col_idx + [1;-1;-1])
            Lfab(i,j) = (n1+1)*(nz+1);
        end
        
        if idx_equals(row_idx, col_idx + [2;-1;-1])
            Lfab(i,j) = -(nz+1)*(nz+2);
        end
        
        if idx_equals(row_idx, col_idx + [0;1;0])
            Sp(i,j) = np + 1;
            Spr(i,j) = np + 1;
        end
        
        if idx_equals(row_idx, col_idx + [-1;1;0])
            Sp(i,j) = -(np + 1);
            Spr(i,j) = (np + 1);
        end
        
        if idx_equals(row_idx, col_idx + [0;0;-1])
            Sp(i,j) = 0.5*(n1 + 1);
            Spr(i,j) = 0.5*(n1 + 1);
        end
        
        if idx_equals(row_idx, col_idx + [1;0;-1])
            Sp(i,j) = 0.5*(nz + 1);
            Spr(i,j) = -0.5*(nz + 1);
        end
        
        if idx_equals(row_idx, col_idx + [0;0;1])
            Sm(i,j) = nm + 1;
            Smr(i,j) = nm + 1;
        end
        
        if idx_equals(row_idx, col_idx + [-1;0;1])
            Sm(i,j) = nm + 1;
            Smr(i,j) = -(nm + 1);
        end
        
        if idx_equals(row_idx, col_idx + [0;-1;0])
            Sm(i,j) = 0.5*(n1 + 1);
            Smr(i,j) = 0.5*(n1 + 1);
        end
        
        if idx_equals(row_idx, col_idx + [1;-1;0])
            Sm(i,j) = -0.5*(nz + 1);
            Smr(i,j) = 0.5*(nz + 1);
        end
        
        row_idx = inc_idx(row_idx, n);
    end
    
    col_idx = inc_idx(col_idx, n);
end

%L = -1i*gamma*gaa/2*Lgaa -1i*gamma*gab/2*Lgab - gamma*fab/2*Lf - w*Lw;
L = -1i*gamma*gaa/2*Lgaa -1i*gamma*gab/2*Lgab - 1i*o/2*(Sp + Sm - Spr - Smr) ...
       -gamma*fab/2*(Sp*Sm + Smr*Spr - 2*Sm*Spr) - gamma*(faa - fab)*Lfaa - w*Lw - gel*Lel;
end

function b = idx_equals(idx1, idx2)
b = min(idx1 == idx2) == 1;
end