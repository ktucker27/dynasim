function [sq2, vars, norms, nps] = get_squeezing(n, es, ess)

sq2 = zeros(size(es,2),1);
vars = zeros(size(es,2),1);
norms = zeros(size(es,2),1);
nps = zeros(3,size(es,2));
for i=1:size(es,2)
    bvn = es(:,i)/norm(es(:,i),2);
    theta = acos(bvn(3));
    if bvn(2) > 0
        phi = acos(bvn(1)/sin(theta));
    else
        phi = 2*pi - acos(bvn(1)/sin(theta));
    end
    
    n1 = [-sin(phi);cos(phi);0];
    n2 = [cos(theta)*cos(phi);cos(theta)*sin(phi);-sin(theta)];
    
    en1n1 = prod_ev(n1, n1, ess(:,:,i));
    en2n2 = prod_ev(n2, n2, ess(:,:,i));
    en1n2 = prod_ev(n1, n2, ess(:,:,i));
    en2n1 = prod_ev(n2, n1, ess(:,:,i));
    
    %G = [en1n1, 0.5*(en1n2 + en2n1);0.5*(en1n2 + en2n1), en2n2];
    %ips(i,:) = eig(G)';
    
    a = en1n1 - en2n2;
    b = en1n2 + en2n1;
    
    phi = 0.5*acos(-a/sqrt(a*a + b*b));
    if b > 0
        phi = pi - phi;
    elseif sqrt(a*a + b*b) == 0
        phi = 0;
    end
    
    np = n1*cos(phi) + n2*sin(phi);
    if(i == 1 || norm(np - nps(:,i-1)) < 0.5)
        nps(:,i) = np;
    else
        nps(:,i) = -np;
    end
    enpnp = prod_ev(np, np, ess(:,:,i));
    
%     ips = zeros(floor(2*pi*100 + 1),2);
%     if true || i == size(es,2)
%         idx = 1;
%         for phi2 = 0:.01:2*pi
%             np2 = n1*cos(phi2) + n2*sin(phi2);
%             ips(idx,1) = prod_ev(np2, np2, ess(:,:,i));
%             idx = idx + 1;
%             if(np2'*es(:,i) > 1.0e-4)
%                 disp('WARNING: Found non-orthogonal direction np2');
%             end
%         end
%         
%         ips(1,2) = enpnp;
%     end
    
    %sq2(i,1) = (1/(2*enpnp))*(en1n1 + en2n2 - sqrt((en1n1 - en2n2)^2 + (en1n2 + en2n1)^2));
    %sq2(i,1) = (2/n)*(en1n1 + en2n2 - sqrt((en1n1 - en2n2)^2 + (en1n2 + en2n1)^2));
    %sq2(i,1) = (n/(2*norm(es(:,i),2)^2))*(en1n1 + en2n2 - sqrt((en1n1 - en2n2)^2 + (en1n2 + en2n1)^2));
    %sq2(i,1) = (n/(2*sum(diag(ess(:,:,i)))))*(en1n1 + en2n2 - sqrt((en1n1 - en2n2)^2 + (en1n2 + en2n1)^2));
    %sq2(i,1) = n*min(real(ips(:,1)))/norm(es(:,i),2)^2;
    %sq2(i,1) = n*min(real(ips(:,1)))/es(1,i)^2;
    sq2(i,1) = n*enpnp/norm(es(:,i),2)^2;
    
    vars(i,1) = enpnp;
    norms(i,1) = (norm(es(:,i),2)^2)/n;
end