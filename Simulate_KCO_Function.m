function [DLCO, KCO] = Simulate_KCO( ...
    RV, TLC, Va_per, LungVols, PulmVol, CapVols, ...
    D0_base, th, th_healthy, tRBC, tRBC_healthy, ...
    th_RBC_barrier, th_scale)

%%FUNCTION [DLCO, KCO] = SIMULATE_KCO(RV, TLC, Va_per, LungVols, PulmVol, CapVols, D0_base, th, th_healthy, tRBC, tRBC_healthy)
%Simulates a simple carbon monoxide gas transfer test using the input parameters
%to generate DLCO and KCO values.
%
%INPUTS:
% RV            - Residual Volume of the lungs (liters)
% TLC           - Total Lung Capacity (liters)
% Va_per        - Percentage of TLC contributing to Alveolar Volume (0-1)
% LungVols      - Relative sizes of different lung regions (array of 4)
% PulmVol       - Total capillary blood volume (liters)
% CapVols      - Relative capillary volumes (array of 4)
% D0_base       - Baseline alveolar-capillary CO diffusion coefficient
% th            - Mean harmonic alveolar-capillary membrane thickness per compartment (μm)
% th_healthy    - Healthy membrane thickness (μm)
% tRBC          - RBC membrane thickness per compartment (nm)
% tRBC_healthy  - Healthy RBC membrane thickness (nm)
% th_RBC_barrier- Permeation analogue (vector of 4)
% th_scale      - Scaling factor for RBC:barrier
%
%OUTPUTS:
% DLCO          - Diffusing capacity of the lungs (mL/min/mmHg)
% KCO           - CO transfer coefficient (mL/min/mmHg/L)

%%Normalise compartment distributions
LungVols = LungVols(:) / sum(LungVols);
CapVols  = CapVols(:) / sum(CapVols);

%%Fixed parameters
BHT = 10;
dt  = 0.001;
C0  = 0.003;

Va = TLC * Va_per;

%%Storage
nSteps = round(BHT / dt);
TotalC_store = zeros(nSteps + 1, 1);

%%Initial concentrations
LungC = C0 * (TLC - RV) * Va_per / Va * ones(4,1);

if abs(TLC - Va) < eps
    Dead_C = 0;
else
    Dead_C = ...
        (C0 * ((TLC - RV) * (1 - Va_per))) / (TLC - Va);
end

TotalC_store(1) = ...
    ((mean(LungC) * Va) + Dead_C * (TLC - Va)) / TLC;

counter = 2;

%%Simulation loop
for t = 0:dt:(BHT-dt)

    for j = 1:4

        %Membrane-thickness pathway
        D0_s = D0_base * ...
            ((((PulmVol * CapVols(j) / 0.065)^(1/2)) / sqrt(2)) ...
            * (Va / TLC)^(1/3)) ...
            * (th_healthy / th(j));

        %Xe-129 RBC:barrier pathway
        %D0_s = D0_base * ...
        %    ((((PulmVol * CapVols(j) / 0.065)^(1/2)) / sqrt(2)) ...
        %    * (Va / TLC)^(1/3)) ...
        %    * th_RBC_barrier(j) * th_scale;

        Pacin = (760 - 47) * LungC(j);
        Pblood = 0;

        Compartment_COAmount = ...
            LungVols(j) * Va * LungC(j);

        newCOAmount = ...
            Compartment_COAmount + ...
            dt * D0_s * (Pblood - Pacin);

        LungC(j) = ...
            newCOAmount / (LungVols(j) * Va);

    end

    TotalC_store(counter) = ...
        ((dot(LungC, LungVols) * Va) + ...
        Dead_C * (TLC - Va)) / TLC;

    counter = counter + 1;
end

%%Calculate DLCO and KCO
Ce = dot(LungC, LungVols);

Cx = C0 * (TLC - RV) * Va_per / Va;

Pb = 760 - 47;

kCO = log(Cx / Ce) / BHT;

DLCO = Va * kCO / Pb;
DLCO = DLCO * 1000 * 60;

KCO = DLCO / (Va * 1.2);

end