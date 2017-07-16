function [Vout, newHint] = quickInterpolate(tVec, zVec, t, hint, round, varargin)
    if nargin < 5
        round = 0;
    end
    if nargin < 4
        hint = 1;
    end
    
    posDir = tVec(1) < tVec(2);
    
    if size(tVec, 1) ~= size(zVec, 1)
        error(message('xDroid:quickInterpolate:dimension not match'));
    end
    
    index = hint;
    
    % looking for index from hint
    % first normalize index
    [index, ~] = rangifyIndex(index, 1, size(tVec, 1)-1);
    isOut = false;
    
    while true
        left = tVec(index);
        right = tVec(index+1);
        
        if posDir
            if left <= t && t < right
                % found index
                break;
            elseif t < left
                index = index - 1;
            elseif t >= right
                index = index + 1;
            end
        else
            if left >= t && t > right
                % found index
                break;
            elseif t > left
                index = index - 1;
            elseif t <= right
                index = index + 1;
            end
        end
        
        [index, isOut] = rangifyIndex(index, 1, size(tVec, 1)-1);
        if isOut
            break;
        end
    end
    newHint = index;
    
    if isOut
        % boundary values, taking as const
        if index == 1
            Vout = zVec(1, :);
        else
            Vout = zVec(end, :);
        end
    else
        % interpolate as usual
        left = tVec(index);
        right = tVec(index+1);
        lambda = (t - left) / (right - left);
        Vout = (1 - lambda) * zVec(index, :) + lambda * zVec(index+1, :);
    end
end

function [index_, isOut] = rangifyIndex(index, l, r)
    if index > r
        index_ = r;
        isOut = true;
    elseif index < l
        index_ = l;
        isOut = true;
    else
        index_ = index;
        isOut = false;
    end
end
