function [minval, finalval, mintime, dur, sstime] = get_min_sq2(t, sq2, sq2ss)

% [~,idx] = max(10*log10(real(sq2)));
%idx = 1;
idx = binary_search(t,.003);
[minval,minidx] = min(10*log10(real(sq2(idx:end))));
idx = minidx + idx - 1;

% [pks, locs] = findpeaks(real(-sq2), t, 'MinPeakProminence', 0.001);
% if size(pks,1) <= 1
%     idx = size(sq2,1);
% else
%     idx = binary_search(t,locs(end,1));
% end
% minval = 10*log10(real(sq2(idx)));
finalval = 10*log10(real(sq2(end,1)));
if minval > finalval
    minval = finalval;
    mintime = t(end,1);
else
    mintime = t(idx);
end

% Calculate duration
% TODO - Assumes increasing function before/after the min with roots
lsq2 = 10*log10(real(sq2));
if idx == size(t,1)
    dur = -1;
else
    idx1 = find_prev_root(lsq2, idx);
    idx2 = find_next_root(lsq2, idx);
    dur = t(idx2) - t(idx1);
end

% Calculate time to steady state
if nargin > 2
    ssidx = find(abs(sq2 - sq2ss)/sq2ss > 0.005, 1, 'last');
    if size(ssidx,1) == 0
        disp('WARNING: Did not reach steady-state');
        sstime = t(end);
    else
        sstime = t(ssidx);
    end
else
    sstime = 0;
end
end

function ridx = find_next_root(v,idx1)
idx2 = size(v,1);
ridx = binary_root(v, idx1, idx2);
end

function ridx = find_prev_root(v,idx2)
[~,idx1] = max(v(1:idx2,1));
ridx = binary_root(v, idx1, idx2);
end

function ridx = binary_root(v, idx1, idx2)
ridx = idx2;
while(abs(v(ridx)) > 5.0e-3 && abs(idx2 - idx1) > 1)
    ridx = floor(0.5*(idx1 + idx2));
    if sign(v(ridx)) == sign(v(idx1))
        idx1 = ridx;
    else
        idx2 = ridx;
    end
end
end