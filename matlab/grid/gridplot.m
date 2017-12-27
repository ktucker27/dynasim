function gridplot(gvec, dvec, A)

figure('Color',[1,1,1])
imagesc(flipud(A))
set(gca, 'YTick', 1:size(dvec,1), 'YTickLabel', flipud(dvec));
set(gca, 'XTick', 1:size(gvec,1), 'XTickLabel', gvec);
colorbar
%xlabel( 'g', 'FontSize', 18 );
%ylabel( '\Delta', 'FontSize', 18 );