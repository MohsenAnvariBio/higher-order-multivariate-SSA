function tplot_figure_13_single(s_k_norm, VR_k, patient_id)
I_dim = length(s_k_norm);
figure('Name', ['Paper Figure 13 - Patient ', num2str(patient_id)], 'Color', 'w', 'Position', [100, 100, 800, 350]);

subplot(1, 2, 1);
plot(1:I_dim, s_k_norm, 'r--', 'LineWidth', 1.5);
xlabel('k', 'FontSize', 12, 'FontAngle', 'italic');
ylabel('Magnitude ||s_k||_2', 'FontSize', 12);
xlim([0, I_dim]);
title('Singular Tube Magnitudes');
grid on;

subplot(1, 2, 2);
plot(1:I_dim, VR_k, 'r--', 'LineWidth', 1.5);
xlabel('k', 'FontSize', 12, 'FontAngle', 'italic');
ylabel('Variance ratio VR(k)', 'FontSize', 12);
xlim([0, I_dim]); ylim([0, 1.05]);
title('Variance Ratio');
grid on;
end
