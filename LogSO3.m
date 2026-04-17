function w = LogSO3(R)
%LOG_SO3 Logarithm map from SO(3) to so(3) (vector form)
%   Input:
%       R - 3x3 rotation matrix
%   Output:
%       w - 3x1 rotation vector

% --- Compute angle safely ---
cos_theta = (trace(R) - 1) / 2;
cos_theta = min(max(cos_theta, -1 + 1e-12), 1 - 1e-12);
theta = acos(cos_theta);

% --- Small-angle case ---
if theta < 1e-8
    w_tilde = 0.5 * (R - R');
    
% --- Near pi case ---
elseif abs(sin(theta)) < 1e-6
    
    RpI = R + eye(3);

    % find most stable column
    col_norms = [norm(RpI(:,1)), norm(RpI(:,2)), norm(RpI(:,3))];
    [~, k] = max(col_norms);

    v = RpI(:,k);
    v = v / norm(v);

    omega = theta * v;

    % skew matrix
    w_tilde = [   0        -omega(3)  omega(2);
               omega(3)      0       -omega(1);
              -omega(2)   omega(1)     0      ];

% --- General case ---
else
    w_tilde = (theta / (2*sin(theta))) * (R - R');
end

% --- vee operator (skew → vector) ---
w = [w_tilde(3,2);
     w_tilde(1,3);
     w_tilde(2,1)];

end