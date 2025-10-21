function plot_filled_arrows_endtips(x, y, col, n_arrows, lenfrac, widthfrac, headfrac)
% Draw filled triangular arrowheads whose tips sit at the END of segments.
% x,y: path vectors. col: [r g b].
% n_arrows   (default: auto): number of heads spaced along path
% lenfrac    (default: 0.05): nominal arrow length as fraction of min axis range
% widthfrac  (default: 0.7) : head width  as fraction of arrow length
% headfrac   (default: 0.8) : head length as fraction of arrow length

if nargin < 4 || isempty(n_arrows), n_arrows = max(1, min(6, floor(numel(x)/6))); end
if nargin < 5 || isempty(lenfrac),   lenfrac   = 0.05; end
if nargin < 6 || isempty(widthfrac), widthfrac = 0.7;  end
if nargin < 7 || isempty(headfrac),  headfrac  = 0.8;  end

dx = diff(x); dy = diff(y);
seg = hypot(dx, dy);
valid = find(seg > 1e-12);
if isempty(valid), return; end

% Use current axes scaling to pick a data-unit arrow size
ax = gca;
Lx = diff(ax.XLim); Ly = diff(ax.YLim);
Lnom = lenfrac * min(Lx, Ly);  % nominal arrow length in data units

% pick segment indices for arrows; the tip will be at index k+1
kpick = valid(round(linspace(1, numel(valid), n_arrows)));
kpick(kpick >= numel(x)) = []; % ensure k+1 exists

% Always add one at the very end
if ~isempty(valid)
    kend = valid(end);
    if kend < numel(x)
        kpick = unique([kpick(:); kend]);  %#ok<AGROW>
    end
end

% Aspect ratio compensation helpers
dar  = daspect(ax);               % [dx dy dz] data units per screen unit
S    = diag(1./dar(1:2));         % to isotropic space
Sinv = diag(dar(1:2));            % back to data space

    for k = kpick(:)'
        % Direction from start of segment to its end (arrow tip at end)
        vx = x(k+1) - x(k);
        vy = y(k+1) - y(k);
        seglen = hypot(vx, vy);
        if seglen <= 1e-12, continue; end
    
        % Use an arrow no longer than the segment (looks cleaner on short steps)
        L = min(Lnom, seglen);
    
        % Unit direction and perpendicular in isotropic space (for nice triangles)
        d_iso = S*[vx; vy]; d_iso = d_iso / norm(d_iso);
        p_iso = [-d_iso(2); d_iso(1)];                  % perpendicular
    
        % Geometry of the head in data coordinates
        tip   = [x(k+1); y(k+1)];
        base  = tip - Sinv*(headfrac*L * d_iso);        % move back along direction
        wvec  = Sinv*((widthfrac*L/2) * p_iso);         % half-width vector
    
        c1 = base + wvec; c2 = base - wvec;
    
        patch([tip(1) c1(1) c2(1)], [tip(2) c1(2) c2(2)], col, ...
              'EdgeColor', col*0.6, 'LineWidth', 0.75, 'FaceAlpha', 1.0);
    end
end

