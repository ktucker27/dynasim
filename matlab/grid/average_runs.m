function [ovec, tf, esf, essf, sq2f, sq2] = average_runs(rootdir, n)

filenames = dir([rootdir, '/*cumulant*_N*.csv']);
pat = '[^_]*cumulant_[^_]*_N[^_]+_D[^_]+_g[^_]+_o([^_])+_f[^_]+_faa[^_]+.csv';

A = zeros(size(filenames,1),2);
for i=1:size(filenames,1)
    t = regexp(filenames(i,1).name, pat, 'tokens');
    A(i,1) = str2double(strrep(t{1}(1), 'p', '.'));
    A(i,2) = i;
    filepaths{i} = [rootdir, '/', filenames(i,1).name];
end

A = sortrows(A,1);
ovec = A(:,1);

total = size(A,1);
for i=1:total
    B = dlmread(filepaths{A(i,2)});
    %[es, ess] = unpack_symm_interp(B, tau);
    [tau, es, ess] = unpack_symm(B);
    
    if i == 1
        tf = tau;
        esf = es;
        essf = ess;
        sq2f = zeros(size(tau,1),1);
        sq2 = zeros(size(tau,1),total);
    else
        if max(abs(tf - tau)) > 1e-10
            disp('ERROR - Received different time vectors');
            return;
        end
        
        esf = esf + es;
        essf = essf + ess;
    end
    
    if nargin > 1
        sq2(:,i) = get_squeezing(n, es, ess);
        %sq2(:,i) = reshape(ess(3,3,:), size(tf));
        %sq2(:,i) = reshape(es(1,:), size(tf));
    end
    
    disp(['o = ', num2str(ovec(i,1))]);
end

esf = esf/total;
essf = essf/total;

if nargin > 1
    sq2f = get_squeezing(n, esf, essf);
    %sq2f = reshape(essf(3,3,:), size(tf));
    %sq2f = reshape(esf(1,:), size(tf));
    
    figure('Color',[1,1,1])
    for i=1:size(sq2,2)
        plot(tf, 10*log10(real(sq2(:,i))), 'Color', [0.8,0.8,0.8])
        %plot(tf, real(sq2(:,i)), 'Color', [0.8,0.8,0.8])
        if i == 1
            hold on
        end
    end
    plot(tf, 10*log10(real(sq2f)), 'r')
    %plot(tf, real(sq2f), 'r')
    xlabel('$\Gamma t$', 'FontSize', 25, 'Interpreter', 'latex')
    ylabel('$\xi^2\mbox{(dB)}$', 'FontSize', 30, 'Interpreter', 'latex')
end