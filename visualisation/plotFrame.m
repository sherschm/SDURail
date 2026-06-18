function plotFrame(T, ls)

    if nargin < 2
        ls = '-';
    end

    origin = T(1:3,4);
    R = T(1:3,1:3);

    scale = 1;

    % X axis (red)
    quiver3(origin(1), origin(2), origin(3), ...
            R(1,1), R(2,1), R(3,1), scale, ...
            'r', 'LineWidth', 2, 'LineStyle', ls);

    %hold on;

    % Y axis (green)
    quiver3(origin(1), origin(2), origin(3), ...
            R(1,2), R(2,2), R(3,2), scale, ...
            'g', 'LineWidth', 2, 'LineStyle', ls);

    % Z axis (blue)
    quiver3(origin(1), origin(2), origin(3), ...
            R(1,3), R(2,3), R(3,3), scale, ...
            'b', 'LineWidth', 2, 'LineStyle', ls);

    %hold off;
end