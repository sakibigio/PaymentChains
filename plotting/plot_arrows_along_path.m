function plot_arrows_along_path(B_path, K_path, col, n_arrows, headscale, lenfrac)
    % Draw small arrows along a parametric 2D path.
    % n_arrows: how many arrows along the path (default auto)
    % headscale: MaxHeadSize for quiver (default 1.5)
    % lenfrac: arrow length as fraction of min axis range (default 0.04)
    
    if nargin < 4 || isempty(n_arrows), n_arrows = max(1, min(6, floor(numel(B_path)/5))); end
    if nargin < 5 || isempty(headscale), headscale = 1.5; end
    if nargin < 6 || isempty(lenfrac), lenfrac = 0.04; end
    
    % Finite differences (segment directions)
    dB = diff(B_path); dK = diff(K_path);
    seglen = hypot(dB, dK);
    valid = find(seglen > 1e-12);
    if isempty(valid), return; end
    
    % Choose arrow locations roughly evenly along non-degenerate segments
    idx = valid(round(linspace(1, numel(valid), n_arrows)));
    
    % Normalize directions and pick a visible arrow length in data units
    ax = gca; Lx = diff(ax.XLim); Ly = diff(ax.YLim);
    L  = lenfrac * min(Lx, Ly);                        % arrow length (data units)
    u  = dB(idx)./seglen(idx) * L;
    v  = dK(idx)./seglen(idx) * L;
    
    % Arrow bases (slightly before the segment end so heads sit on the curve)
    x0 = B_path(idx);
    y0 = K_path(idx);
    
    % Draw arrows (no autoscale => exact lengths)
    quiver(x0, y0, u, v, 0, ...
        'Color', col, 'LineWidth', 1.5, 'MaxHeadSize', headscale);
    
    % Also put a prominent arrow at the very end, if there is motion
    iend = valid(end);
    ue = dB(iend)/seglen(iend)*L; ve = dK(iend)/seglen(iend)*L;
   % quiver(B_path(iend), K_path(iend), ue, ve, 0, ...
     %   'Color', col, 'LineWidth', 1.5, 'MaxHeadSize', 2.2);
end
