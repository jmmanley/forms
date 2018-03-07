function ts = sil_sample(res, ctrl_pts, error)

% SIL_SAMPLE  Return parameter values to evenly sample silhouette given by
%   ctrl_pts to resolution res.

if nargin < 3
    error = 1e-6;
end

lengths = [];
splits  = [];

for i = 1:size(ctrl_pts, 3)
    seg_ctrl = ctrl_pts(:, :, i);
    
    [seg_lengths seg_splits] = sil_flattenbezier(seg_ctrl, error);
    
    lengths = [ lengths seg_lengths            ];              %#ok<AGROW>
    splits  = [ splits  (i - 1) + seg_splits i ];              %#ok<AGROW>
end

total_length = sum(lengths);

cum_lengths = cumsum(lengths);
length_params = linspace(0, total_length * (res - 1) / res, res);
ts = zeros(1, res);

for i = 1:res
    end_seg = find(cum_lengths > length_params(i), 1);
    
    end_length = cum_lengths(end_seg);
    end_t = splits(end_seg);
    
    if end_seg == 1
        start_length = 0;
        start_t = 0;
    else
        start_length = cum_lengths(end_seg - 1);
        start_t = splits(end_seg - 1);
    end
    
    ts(i) = (length_params(i) - start_length) * (end_t - start_t) / ...
            (end_length - start_length) + start_t;
end

end