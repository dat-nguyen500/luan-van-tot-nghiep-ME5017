clc; clear; close all;

%% ====== THÔNG SỐ CƠ CẤU ======
% Mục tiêu script:
%   - Mô phỏng cơ cấu phẳng với khâu dẫn OA quay đều.
%   - Tại mỗi thời điểm, giải hệ hình học để tìm điểm B (giao 2 đường tròn).
%   - Từ A, B suy ra điểm M theo quan hệ M = 2B - A.
%   - Vẽ animation cơ cấu và đồng thời vẽ (dạng điểm) quan hệ tuyến tính gần đúng giữa c_x và theta.
%   - Sau khi mô phỏng xong: tính metric % cho c_y và v_x, và fit mô hình tuyến tính: c_x = a*theta + b.
%
% Lưu ý đơn vị:
%   - Tọa độ (mm), góc (rad), thời gian (s), vận tốc (mm/s)
O = [0, 0];
C = [54, 0];

% Chiều dài các khâu/ bán kính ràng buộc hình học
r1 = 27;   % OA
r2 = 67;   % AB (ràng buộc: |B-A| = r2)
r3 = 67;   % BC (ràng buộc: |B-C| = r3)

% Vận tốc góc của khâu dẫn (OA) - thay đổi omega1 ở đây
omega1 = 5;              % rad/s

% Miền làm việc của góc theta (rad)
% Ở đây chọn quét từ theta_min đến theta_max.
theta0 = 7*pi/9;           % góc ban đầu tại t = 0
theta_min = 7*pi/9;        % góc nhỏ nhất
theta_max = 11*pi/9;        % góc lớn nhất

% Thời gian kết thúc tương ứng với quét góc [theta_min, theta_max]
t_end = (theta_max-theta_min)/omega1;
% Chọn số mẫu thời gian (càng lớn animation càng mượt nhưng chậm hơn)
t = linspace(0, t_end, 200);

%% ====== GIẢI HỆ HÌNH HỌC (SYMBOLIC) + KHAI BÁO BỘ NHỚ ======
% Dùng symbolic để giải giao điểm 2 đường tròn (tìm tọa độ B).
% x,y là biến symbolic.
syms x y real

% B_prev: dùng để chọn nghiệm liên tục của B (tránh nhảy giữa 2 nghiệm giao nhau)
B_prev = [];

% Pre-allocate để tăng tốc và tránh nối mảng trong vòng lặp
N = length(t);
M_trace = nan(N, 2);
theta_hist = nan(N, 1);
cx_hist = nan(N, 1);
cy_hist = nan(N, 1);
vx_hist = nan(N, 1);

%% ====== FIGURE/AXES ======
% Gộp mô phỏng + đồ thị c_x-theta vào 1 cửa sổ duy nhất
figMain = figure;
tiled = tiledlayout(figMain, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axSim = nexttile(tiled, 1);
% Axes mô phỏng cơ cấu
axis(axSim, 'equal'); grid(axSim, 'on'); hold(axSim, 'on')
xlim(axSim, [-40 100])
ylim(axSim, [-20 160])
xlabel(axSim, 'x (mm)')
ylabel(axSim, 'y (mm)')

axThetaX = nexttile(tiled, 2);
grid(axThetaX, 'on'); hold(axThetaX, 'on')
% Đồ thị quan hệ c_x - theta: vẽ dạng điểm (không nối đoạn thẳng)
hThetaX = plot(axThetaX, nan, nan, 'm.', 'LineStyle', 'none', 'MarkerSize', 10);
xlabel(axThetaX, '\theta_1 (rad)')
ylabel(axThetaX, 'c_x (mm)')
title(axThetaX, 'Đồ thị quan hệ c_x - \theta_1')

%% ====== QUỸ ĐẠO KHÂU DẪN OA (VẼ TRƯỚC ANIMATION) ======
% Yêu cầu: vẽ trước quỹ đạo OA rồi mới chạy animation.
% - 2 đoạn nét đứt: OA ở vị trí đầu và vị trí cuối
% - 1 cung tròn nét đứt: quỹ đạo điểm A (do A quay quanh O bán kính r1)
theta_arc = linspace(theta_min, theta_max, 200);
A_arc = [r1*cos(theta_arc(:)), r1*sin(theta_arc(:))];
A_start = [r1*cos(theta_min), r1*sin(theta_min)];
A_end = [r1*cos(theta_max), r1*sin(theta_max)];

plot(axSim, [O(1) A_start(1)], [O(2) A_start(2)], 'k--', 'LineWidth', 1)
plot(axSim, [O(1) A_end(1)], [O(2) A_end(2)], 'k--', 'LineWidth', 1)
plot(axSim, A_arc(:,1), A_arc(:,2), 'k--', 'LineWidth', 1)
drawnow

%% ====== MÔ PHỎNG ======
% hDyn: tập handle các đối tượng "động" (vẽ lại mỗi frame).
% Ta xóa hDyn mỗi vòng để giữ lại các đối tượng nền (quỹ đạo OA nét đứt).
hDyn = gobjects(0);
for k = 1:length(t)

    % Góc khâu dẫn OA tại thời điểm t(k)
    theta = theta0 + omega1 * t(k);

    % Tọa độ điểm A trên đường tròn tâm O bán kính r1
    A = [r1*cos(theta), r1*sin(theta)];

    % Tìm điểm B là giao của 2 đường tròn:
    %   (1) tâm A, bán kính r2
    %   (2) tâm C, bán kính r3
    % => giải hệ 2 phương trình để ra (x,y) = B
    eqn1 = (x-A(1))^2 + (y-A(2))^2 == r2^2;
    eqn2 = (x-C(1))^2 + (y-C(2))^2 == r3^2;

    % Lấy nghiệm thực (thường có 2 nghiệm đối xứng qua đường AC)
    sol = solve([eqn1, eqn2], [x y], 'Real', true);
    B_candidates = double([sol.x, sol.y]);

    % Chọn nghiệm B liên tục theo thời gian để tránh nhảy nghiệm:
    % - Frame đầu: chọn nghiệm phía trên (y > 0)
    % - Các frame sau: chọn nghiệm gần B_prev nhất
    if k == 1
        idx = B_candidates(:,2) > 0;
        B = B_candidates(idx,:);
    else
        d = vecnorm(B_candidates - B_prev, 2, 2);
        [~, idx] = min(d);
        B = B_candidates(idx,:);
    end

    B_prev = B;
    % Từ A, B suy ra điểm M theo quan hệ hình học của cơ cấu
    % (đang dùng theo yêu cầu: M = 2B - A)
    M = 2*B - A;

    % Lưu lịch sử để:
    % - vẽ quỹ đạo M
    % - vẽ quan hệ c_x - theta
    % - tính metric % sau mô phỏng
    theta_hist(k) = theta;
    cx_hist(k) = M(1);
    cy_hist(k) = M(2);
    M_trace(k,:) = M;

    % Vận tốc theo phương x của điểm M (không quan tâm chiều => lấy trị tuyệt đối)
    % Dùng sai phân hữu hạn bậc 1: v_x(k) = | (c_x(k) - c_x(k-1)) / dt |
    if k == 1
        vx_hist(k) = NaN;
    else
        dt = t(k) - t(k-1);
        % Lấy tốc độ theo phương x (không quan tâm chiều)
        vx_hist(k) = abs((cx_hist(k) - cx_hist(k-1)) / dt);
    end

    %% Vẽ (giữ lại quỹ đạo OA đã vẽ sẵn)
    % Lưu ý: Không dùng cla() vì sẽ xóa cả nền quỹ đạo OA nét đứt.
    % Thay vào đó, xóa các đối tượng động của frame trước bằng handle hDyn.
    if ~isempty(hDyn)
        delete(hDyn(isgraphics(hDyn)))
    end

    % Vẽ các thanh liên kết
    h1 = plot(axSim, [O(1) A(1)], [O(2) A(2)], 'r', 'LineWidth', 2); % khâu dẫn OA
    h2 = plot(axSim, [A(1) M(1)], [A(2) M(2)], 'g', 'LineWidth', 2);
    h3 = plot(axSim, [B(1) C(1)], [B(2) C(2)], 'b', 'LineWidth', 2);

    % Vẽ các điểm nút
    hO = plot(axSim, O(1), O(2), 'ko','MarkerFaceColor','k');
    hA = plot(axSim, A(1), A(2), 'ro','MarkerFaceColor','r');
    hB = plot(axSim, B(1), B(2), 'go','MarkerFaceColor','g');
    hC = plot(axSim, C(1), C(2), 'bo','MarkerFaceColor','b');
    hM = plot(axSim, M(1), M(2), 'mo','MarkerFaceColor','m');       % điểm M

    % Gán tên các điểm để quan sát trực quan (vị trí chữ đã được chỉnh theo yêu cầu)
    labelOffset = [2, 2];
    tO = text(axSim, O(1)+labelOffset(1), O(2)+labelOffset(2), 'O', 'FontWeight', 'bold');
    tB = text(axSim, B(1)+labelOffset(1), B(2)+labelOffset(2), 'B', 'FontWeight', 'bold');
    tC = text(axSim, C(1)+labelOffset(1), C(2)+labelOffset(2), 'C', 'FontWeight', 'bold');

    % A nằm bên trái điểm A
    tA = text(axSim, A(1)-4, A(2), 'A', 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    % M nằm phía trên điểm M
    tM = text(axSim, M(1), M(2)+4, 'M', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    % Quỹ đạo điểm M (đã chạy tới frame k)
    hTrace = plot(axSim, M_trace(1:k,1), M_trace(1:k,2), 'm', 'LineWidth', 1.5);

    hDyn = [h1, h2, h3, hO, hA, hB, hC, hM, tO, tA, tB, tC, tM, hTrace];

    % Tiêu đề mô phỏng (đơn giản, chỉ hiển thị omega)
    title(axSim, sprintf('Mô phỏng với \\omega_1 = %.2f rad/s', omega1))

    % Cập nhật đồ thị quan hệ (theta, c_x) theo thời gian
    set(hThetaX, 'XData', theta_hist(1:k), 'YData', cx_hist(1:k));
    drawnow limitrate
end

%% ====== HỒI QUY TUYẾN TÍNH: c_x = a*theta + b ======
% Mục tiêu: tìm mô hình tuyến tính gần đúng giữa c_x và theta theo dạng:
%   c_x = a*theta + b
% Dùng least-squares thông qua polyfit bậc 1.
fitMask = isfinite(theta_hist) & isfinite(cx_hist);
theta_fit = theta_hist(fitMask);
cx_fit = cx_hist(fitMask);

if numel(theta_fit) >= 2
    % p = [a, b] sao cho c_x ≈ a*theta + b
    p = polyfit(theta_fit, cx_fit, 1);
    a = p(1);
    b = p(2);

    % Dự đoán và sai số mô hình
    cx_pred = a*theta_fit + b;
    resid = cx_fit - cx_pred;

    % Các metric khớp mô hình (đơn vị mm hoặc vô hướng)
    rmse_cx = sqrt(mean(resid.^2));
    mae_cx = mean(abs(resid));
    max_abs_err_cx = max(abs(resid));

    % Hệ số xác định R^2
    ss_res = sum(resid.^2);
    ss_tot = sum((cx_fit - mean(cx_fit)).^2);
    if ss_tot == 0
        r2 = NaN;
    else
        r2 = 1 - ss_res/ss_tot;
    end

    % Sai số phần trăm cực đại (chuẩn hóa theo mean(|c_x|) để không phụ thuộc dấu)
    denom = mean(abs(cx_fit));
    if denom == 0
        max_err_pct = NaN;
    else
        max_err_pct = max_abs_err_cx * 100 / denom;
    end
else
    a = NaN;
    b = NaN;
    rmse_cx = NaN;
    mae_cx = NaN;
    max_abs_err_cx = NaN;
    r2 = NaN;
    max_err_pct = NaN;
end

%% ====== METRIC (Sai số phần trăm so với giá trị trung bình) ======
% Theo yêu cầu:
%   err_cy(%) = max(abs(cy - mean(cy))) * 100 / abs(mean(cy))
%   err_vx(%) = max(abs(vx - mean(vx))) * 100 / abs(mean(vx))
% Trong đó:
%   - cy là tọa độ y của điểm M (đặt tên cy cho báo cáo)
%   - vx là tốc độ theo x của điểm M (đã lấy trị tuyệt đối nên luôn dương)
cy_mean = mean(cy_hist, 'omitnan');
cy_dev_max = max(abs(cy_hist - cy_mean), [], 'omitnan');
cy_denom = abs(cy_mean);
if cy_denom == 0
    cy_err_pct = NaN;
else
    cy_err_pct = cy_dev_max * 100 / cy_denom;
end

vx_valid = vx_hist(isfinite(vx_hist));
vx_mean = mean(vx_valid);
vx_dev_max = max(abs(vx_valid - vx_mean));
vx_denom = abs(vx_mean);
if vx_denom == 0
    vx_err_pct = NaN;
else
    vx_err_pct = vx_dev_max * 100 / vx_denom;
end

fprintf('\n===== ĐÁNH GIÁ (%%) =====\n');
fprintf('mean(cy) = %.6f mm\n', cy_mean);
fprintf('err_cy(%%) = max(abs(cy - mean(cy))) * 100 / abs(mean(cy)) = %.6f %%\n', cy_err_pct);
fprintf('mean(vx) = %.6f mm/s\n', vx_mean);
fprintf('err_vx(%%) = max(abs(vx - mean(vx))) * 100 / abs(mean(vx)) = %.6f %%\n', vx_err_pct);

fprintf('\n===== MÔ HÌNH TUYẾN TÍNH =====\n');
fprintf('c_x = a*theta + b\n');
fprintf('a = %.6f (mm/rad)\n', a);
fprintf('b = %.6f (mm)\n', b);
fprintf('R^2 = %.6f\n', r2);
fprintf('RMSE(c_x) = %.6f mm\n', rmse_cx);
fprintf('MAE(c_x) = %.6f mm\n', mae_cx);
fprintf('max|err(c_x)| = %.6f mm\n', max_abs_err_cx);
fprintf('max|err(c_x)| (%% of mean|c_x|) = %.6f %%\n', max_err_pct);

fprintf('\n===== HÀNH TRÌNH c_x =====\n');
fprintf('c_x(start) = %.6f mm\n', cx_hist(1));
fprintf('c_x(end)   = %.6f mm\n', cx_hist(end));
fprintf('Δc_x = c_x(end) - c_x(start) = %.6f mm\n', cx_hist(end) - cx_hist(1));
fprintf('|Δc_x| = %.6f mm\n', abs(cx_hist(end) - cx_hist(1)));
