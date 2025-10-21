% === Auto-sync figures to Overleaf ===
src = 'Matlab Figures/';
dst = '/Users/sakiclaudia/Dropbox/Apps/Overleaf/Chained Payments/Matlab Figures/';

if ~exist(dst, 'dir')
    mkdir(dst);
end

extensions = {'*.pdf', '*.eps'};

for k = 1:length(extensions)
    files = dir(fullfile(src, extensions{k}));
    for i = 1:length(files)
        copyfile(fullfile(src, files(i).name), fullfile(dst, files(i).name));
    end
end

disp('Figures copied to Overleaf.');
