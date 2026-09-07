function [DLCO, KCO] = Simulate_KCO_qCT( ...
    RV, TLC, Va_per, LungVols, PulmVol, CapVols, ...
    D0_base, th, th_healthy, tRBC, tRBC_healthy, ...
    th_RBC_barrier, th_scale, SA_scale)

%Simulate_KCO_qCT
%Simulates a 10-s single-breath CO gas-transfer manoeuvre using
%qCT-derived regional structural inputs.
%
%INPUTS:
% RV             - Residual volume [L]
% TLC            - Total lung capacity [L]
% Va_per         - Fraction of TLC contributing to alveolar volume [0-1]
% LungVols       - Regional fractional distribution of alveolar volume
% PulmVol        - Total pulmonary capillary blood volume [L]
% CapVols        - Regional capillary blood-volume distribution
% D0_base        - Subject-specific baseline diffusion coefficient
% th             - ACM thickness per compartment [micrometres]
% th_healthy     - Healthy reference ACM thickness [micrometres]
% tRBC           - RBC membrane thickness per compartment [nm]
% tRBC_healthy   - Healthy reference RBC membrane thickness [nm]
% th_RBC_barrier - RBC:barrier scaling term
% th_scale       - Additional RBC:barrier scaling factor
% SA_scale       - Relative functional gas-exchange surface-area scaling
%
%OUTPUTS:
% DLCO           - mL/min/mmHg
% KCO            - mL/min/mmHg/L

%%Normalise regional distributions

LungVols = LungVols(:) / sum(LungVols);
CapVols  = CapVols(:)  / sum(CapVols);

SA_scale = SA_scale(:);

if numel(SA_scale) == 1
    SA_scale = repmat(SA_scale, 4, 1);
end

%Fixed simulation parameters

BHT = 10;          %Breath-hold time [s]
dt  = 0.001;       %Timestep [s]
C0  = 0.003;       %Initial inspired CO fraction

Va = TLC * Va_per;

%%Initialise storage

nSteps = round(BHT / dt);
TotalC_store = zeros(nSteps + 1, 1);

%%Initial alveolar concentrations

LungC = ...
    C0 * (TLC - RV) * Va_per / Va * ones(4,1);

if abs(TLC - Va) < eps

    Dead_C = 0;

else

    Dead_C = ...
        (C0 * ((TLC - RV) * (1 - Va_per))) / ...
        (TLC - Va);

end

TotalC_store(1) = ...
    ((mean(LungC) * Va) + ...
    Dead_C * (TLC - Va)) / TLC;

counter = 2;

%%Simulation loop

for t = 0:dt:(BHT-dt)

    for j = 1:4

        %qCT membrane-thickness + surface-area pathway
        D0_s = D0_base * ...
            (((PulmVol * CapVols(j) / 0.065)^(1/2) / sqrt(2)) ...
            * (Va / TLC)^(1/3)) ...
            * (th_healthy / th(j)) ...
            * SA_scale(j);

        %Alveolar CO partial pressure [mmHg]
        Pacin = (760 - 47) * LungC(j);

        %Pulmonary capillary blood approximated as an infinite CO sink
        Pblood = 0;

        %Current compartment CO amount
        Compartment_COAmount = ...
            LungVols(j) * Va * LungC(j);

        %CO uptake over timestep
        newCOAmount = ...
            Compartment_COAmount + ...
            dt * D0_s * (Pblood - Pacin);

        %Updated alveolar concentration
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

%Convert L/s/mmHg to mL/min/mmHg
DLCO = DLCO * 1000 * 60;

%Retain original study implementation
KCO = DLCO / (Va * 1.2);

end