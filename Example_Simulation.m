%example_simulation
RV = 1.5;
TLC = 6.0;
Va_per = 0.85;

LungVols = [0.23 0.237 0.264 0.269];
PulmVol = 0.065;
CapVols = [0.25 0.25 0.25 0.25];

D0_base = 0.001;

th = [1 1 1 1];
th_healthy = 1;

tRBC = [0.15 0.15 0.15 0.15];
tRBC_healthy = 0.15;

th_RBC_barrier = [1 1 1 1];
th_scale = 1;

[DLCO, KCO] = Simulate_KCO( ...
    RV, TLC, Va_per, LungVols, PulmVol, CapVols, ...
    D0_base, th, th_healthy, tRBC, tRBC_healthy, ...
    th_RBC_barrier, th_scale);

fprintf('Simulated DLCO: %.2f mL/min/mmHg\n', DLCO);
fprintf('Simulated KCO: %.2f mL/min/mmHg/L\n', KCO);