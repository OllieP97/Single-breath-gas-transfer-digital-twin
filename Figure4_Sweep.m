%%FIGURE 4 SIMULATIONS
%Reproduces the mechanistic simulations shown in Figure 4.
%
%Requires:
%  Simulate_KCO.m

clear;
clc;
close all;


%BASE SIMULATION

RV = 1.5;                  %Residual volume [L]
TLC = 6.0;                 %Total lung capacity [L]
Va_per = 0.85;             %Fraction of TLC contributing to alveolar volume

%Four compartments: LUL, LLL, RUL+RML, RLL
LungVols = [0.23 0.237 0.264 0.269];

PulmVol = 0.065;           %Pulmonary capillary blood volume [L]
CapVols = [0.25 0.25 0.25 0.25];

D0_base = 0.001;

th = [1.0 1.0 1.0 1.0];   %ACM thickness [micrometres]
th_healthy = 1.0;

tRBC = [0.15 0.15 0.15 0.15];
tRBC_healthy = 0.15;

th_RBC_barrier = [1 1 1 1];
th_scale = 1;

%%Run baseline simulation

[DLCO, KCO] = Simulate_KCO( ...
    RV, TLC, Va_per, LungVols, PulmVol, CapVols, ...
    D0_base, th, th_healthy, tRBC, tRBC_healthy, ...
    th_RBC_barrier, th_scale);

fprintf('Baseline simulated DLCO: %.2f mL/min/mmHg\n', DLCO);
fprintf('Baseline simulated KCO: %.2f mL/min/mmHg/L\n', KCO);




%VENTILATION HETEROGENEITY: LOBAR COEFFICIENT OF VARIATION

rng(1);                           % Reproducible random sampling

num_samples = 300;
Va_vals = linspace(0.7, 0.95, 5);

DLCO_mat = NaN(num_samples, numel(Va_vals));
KCO_mat  = NaN(num_samples, numel(Va_vals));

CV_vals = NaN(num_samples, 1);
LungVols_all = NaN(num_samples, 4);

%Physiological supine regional ventilation ranges
%Compartments: [LUL, LLL, RUL+RML, RLL]

vlow  = [0.15 0.15 0.15 0.15];
vhigh = [0.30 0.30 0.35 0.35];

%Generate ventilation distributions

for i = 1:num_samples

    %Sample compartmental weights from physiological ranges
    x = vlow + (vhigh - vlow) .* rand(1,4);

    %Normalise regional ventilation distribution to unity
    x = x / sum(x);

    LungVols_all(i,:) = x;

    %Coefficient of variation across four compartments
    CV_vals(i) = std(x) / mean(x);

end

%%Run Va/TLC x heterogeneity sweep

for i = 1:num_samples

    LungVols_i = LungVols_all(i,:);

    for j = 1:numel(Va_vals)

        Va_per_j = Va_vals(j);

        [DLCO_ij, KCO_ij] = Simulate_KCO( ...
            RV, TLC, Va_per_j, LungVols_i, ...
            PulmVol, CapVols, D0_base, ...
            th, th_healthy, ...
            tRBC, tRBC_healthy, ...
            th_RBC_barrier, th_scale);

        DLCO_mat(i,j) = DLCO_ij;
        KCO_mat(i,j)  = KCO_ij;

    end
end

%%Sort according to ventilation CV

[CV_sorted, idx] = sort(CV_vals);

DLCO_CV = DLCO_mat(idx,:);
KCO_CV  = KCO_mat(idx,:);

%%Displayed CV range

cv_min = max(0, prctile(CV_vals,1) - 0.005);
cv_max = prctile(CV_vals,99) + 0.005;

[CV_grid, VA_grid] = meshgrid(CV_sorted, Va_vals);

%%KCO surface

figure;

surf( ...
    CV_sorted, ...
    Va_vals, ...
    KCO_CV', ...
    'EdgeColor','none');

xlabel('Lobar ventilation CV');
ylabel('V_A/TLC');
zlabel('KCO (mL/min/mmHg/L)');

title('KCO vs ventilation heterogeneity and V_A/TLC');

xlim([cv_min cv_max]);

shading interp;
colorbar;
hold on;

contour3( ...
    CV_grid, ...
    VA_grid, ...
    KCO_CV', ...
    20, ...
    'k');

view(3);


%%DLCO surface

figure;

surf( ...
    CV_sorted, ...
    Va_vals, ...
    DLCO_CV', ...
    'EdgeColor','none');

xlabel('Lobar ventilation CV');
ylabel('V_A/TLC');
zlabel('DLCO (mL/min/mmHg)');

title('DLCO vs ventilation heterogeneity and V_A/TLC');

xlim([cv_min cv_max]);

shading interp;
colorbar;
hold on;

contour3( ...
    CV_grid, ...
    VA_grid, ...
    DLCO_CV', ...
    20, ...
    'k');

view(3);



%UPPER-TO-LOWER VENTILATION REDISTRIBUTION
Va_vals_UL = linspace(0.5, 0.9, 20);

%Avoid exactly 0 and 1 because these create zero-volume compartments
upper_fracs = linspace(0.01, 0.99, 20);

DLCO_UL = zeros(length(upper_fracs), length(Va_vals_UL));
KCO_UL  = zeros(length(upper_fracs), length(Va_vals_UL));

%%Run redistribution sweep

for i = 1:length(upper_fracs)

    upper = upper_fracs(i);
    lower = 1 - upper;

    %Regional fractional ventilation:
    %[LUL, LLL, RUL+RML, RLL]
    %
    %Total upper ventilation is divided equally between the left and right upper compartments; lower ventilation is divided equally between the left and right lower compartments.

    LungVols_UL = [ ...
        upper/2, ...
        lower/2, ...
        upper/2, ...
        lower/2];

    for j = 1:length(Va_vals_UL)

        Va_per_j = Va_vals_UL(j);

        [DLCO_UL(i,j), KCO_UL(i,j)] = Simulate_KCO( ...
            RV, TLC, Va_per_j, LungVols_UL, ...
            PulmVol, CapVols, D0_base, ...
            th, th_healthy, ...
            tRBC, tRBC_healthy, ...
            th_RBC_barrier, th_scale);

    end
end

%%Plotting grid
%Convert upper fraction from 0-1 to percentage for display.

[UF, VA_UL] = meshgrid( ...
    upper_fracs * 100, ...
    Va_vals_UL);


%%DLCO: upper-to-lower redistribution

figure;

surf( ...
    UF, ...
    VA_UL, ...
    DLCO_UL', ...
    'EdgeColor','none');

xlabel('Ventilation to upper lung (%)');
ylabel('V_A/TLC');
zlabel('DLCO (mL/min/mmHg)');

title('DLCO vs upper-to-lower ventilation redistribution');

shading interp;
colorbar;
view(3);


%%KCO: upper-to-lower redistribution

figure;

surf( ...
    UF, ...
    VA_UL, ...
    KCO_UL', ...
    'EdgeColor','none');

xlabel('Ventilation to upper lung (%)');
ylabel('V_A/TLC');
zlabel('KCO (mL/min/mmHg/L)');

title('KCO vs upper-to-lower ventilation redistribution');

shading interp;
colorbar;
view(3);





%ACM Thickness Plots

%Field plot
thickness_vals = linspace(0.5, 3, 20);        %membrane thickness range (μm)
pulmvol_vals = linspace(0.03, 0.1, 20);       %pulmonary blood volume range (liters)

DLCO_mat = zeros(length(thickness_vals), length(pulmvol_vals));
KCO_mat = zeros(length(thickness_vals), length(pulmvol_vals));

RV = 1.5;
TLC = 6.0;
Va_per = 0.85;
LungVols = [0.23 0.237 0.264 0.269];
CapVols = [0.25 0.25 0.25 0.25];
D0_base = 0.001;
th_healthy = 1.0;
tRBC = [0.15 0.15 0.15 0.15];
tRBC_healthy = 0.15;
th_RBC_barrier = [1 1 1 1];
th_scale = 1;

for i = 1:length(thickness_vals)
    for j = 1:length(pulmvol_vals)
        th = thickness_vals(i) * ones(1,4);        %set all compartments same
        PulmVol = pulmvol_vals(j);
        [DLCO, KCO] = Simulate_KCO(RV, TLC, Va_per, LungVols, PulmVol, CapVols, ...
                                   D0_base, th, th_healthy, tRBC, tRBC_healthy, ...
                                   th_RBC_barrier, th_scale);
        DLCO_mat(i,j) = DLCO;
        KCO_mat(i,j) = KCO;
    end
end

[TH, PV] = meshgrid(thickness_vals, pulmvol_vals);

%Plot DLCO
figure;
surf(TH, PV, DLCO_mat')
xlabel('Membrane Thickness (μm)')
ylabel('Pulmonary Blood Volume (L)')
zlabel('DLCO (mL/min/mmHg)')
title('DLCO variation with thickness and blood volume')
colorbar
shading interp

%Plot KCO
figure;
surf(TH, PV, KCO_mat')
xlabel('Membrane Thickness (μm)')
ylabel('Pulmonary Blood Volume (L)')
zlabel('KCO (mL/min/mmHg/L)')
title('KCO variation with thickness and blood volume')
colorbar
shading interp