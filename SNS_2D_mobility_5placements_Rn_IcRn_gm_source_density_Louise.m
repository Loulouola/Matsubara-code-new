%% SNS_2D_mobility_5placements_Rn_IcRn_gm.m
% Model overview:

% Wide 2D SNS junction: k_y continuum + exact parabolic BdG Matsubara model.
%
% Mobility is calculated locally from the spatially varying sheet density:
%   mu_mob(x) = mu_ref * [n_s(x)/n_ref]^alpha_mob
%   lambda_mfp(x) = hbar*mu_mob(x)*sqrt(2*pi*n_s(x))/e
%                = lambda_ref * [n_s(x)/n_ref]^(alpha_mob+1/2)
%
% Incoherent distributed resistance is replaced by ONE effective scatterer:
%   rho = integral_0^L dx/lambda(x) % dimensionless measure of accumulated scattering
%   tau_tot = 1/(1+rho) % effective transmssion
%
% The SAME tau_tot applies to all k_y channels at a chosen mu_N. The coherent
% scatterer is inserted at 5 fixed, randomly selected positions in [0,L].
% Its location affects quantum interference but NOT the incoherent tau_tot.
%
% The primary result is mean_j{max_phi |I_j(phi)|}, NOT
% max_phi |mean_j I_j(phi)|. Both are available for comparison.
%
% Models: smooth band profile, abrupt band profile, ideal ballistic Andreev.
% Includes normal-state G_N, R_N, Ic*R_N and dIc/dVg via Cox.
% Source chemical potential is determined from source 2D sheet density
% and source effective mass (one occupied subband, spin degeneracy 2).
% All spatial integration includes the full normal length
% Lch = Lg + 2*Lungated, including both underlaps.
%
% This is a phenomenological (simple effective rule) model. The fixed-tau point scatterer is NOT a
% literal delta-potential and its analytic continuation to imaginary energies
% is an assumption. This first step uses the LOCAL 2D Thomas-Fermi density,
% justified here for mu_N >= 2*Delta, rather than quantum charge corrections.
%
% MATLAB R2020a or newer (uses tiledlayout and script-local functions).

% No global clear/close here so this packaged runner is safe to call from a wrapper.

%% Physical constants
% p is a parameter structure: a bundle of simulation settings and physical constants. 
p.e    = 1.602176634e-19; % elemental charge
p.hbar = 1.054571817e-34; % Reduced Planck's constant
p.kB   = 1.380649e-23; % Boltzmann's constant
p.m0   = 9.1093837139e-31; % electron mass
p.eV   = p.e; % electron volt
p.h    = 2*pi*p.hbar; % Planck's constant

%% Geometry: same Lg / underlap / electrostatic-ramp construction as Fig. 1
% The Fig. 1 illustration values are used as the single-geometry default.
% Set Lg_nm to another valid value to study a different gate length.
p.Lg_nm = 40;                       % Physical gate-foot length [nm]
lgOverride=getenv('LOUISE_LG_NM');
if ~isempty(lgOverride), p.Lg_nm=str2double(lgOverride); end
p.Lungated_nm = 20;                 % Ungated underlap on EACH side [nm]
p.tox_nm = 12;                      % Gate dielectric thickness [nm]
toxOverride=getenv('LOUISE_TOX_NM');
if ~isempty(toxOverride)
    p.tox_nm=str2double(toxOverride);
    assert(isfinite(p.tox_nm) && p.tox_nm>0, ...
        'LOUISE_TOX_NM must be a positive finite value in nm.');
end
p.tch_nm = 12;                      % Effective full channel thickness [nm]
p.alphaGeo = 1;                     % Electrostatic boundary-condition factor
p.epsilonChannelRelative = 16.4;    % Channel relative permitttivity
p.epsilonOxRelative = 9.0;          % Oxide relative permittivity
p.epsilon0 = 8.8541878128e-12;      % Permittivity of vacuum
p.lambda_nm = sqrt(p.alphaGeo*p.epsilonChannelRelative / ...
    p.epsilonOxRelative*p.tch_nm*p.tox_nm); % "electrostatic transition length
p.Lch_nm = p.Lg_nm + 2*p.Lungated_nm; % Channel length
p.Ltrans_nm = p.Lungated_nm + p.lambda_nm; % Transition length for chemical potential
p.Lplateau_nm = p.Lg_nm - 2*p.lambda_nm; % Platau length with constant chemical potential
p.L = p.Lch_nm*1e-9;                % Entire normal region [m]
p.Ltrans = p.Ltrans_nm*1e-9;        % Smooth transition on EACH side [m]

%% Gate capacitance and voltage conversion
% Use the same oxide permittivity and thickness as the geometry above.
p.Cox = p.epsilon0*p.epsilonOxRelative/(p.tox_nm*1e-9);
p.Vfb = 0;                 % Adjustable effective flat-band offset [V]

%% Material / device parameters
p.mN = 0.023*p.m0;                  % Effective mass in normal region
p.mS = 0.023*p.m0;                % Source/S-region effective mass [kg]
p.nSource_cm2 = 2.0e12;        % Source sheet density [cm^-2] INPUT
p.nSource = p.nSource_cm2*1e4;  % Source sheet density [m^-2]
p.gValleySource = 1;            % Source valley degeneracy
% Single occupied 2D subband; spin degeneracy g_s=2:
% nS = g_v*mS*muS/(pi*hbar^2) at T=0.
p.muS = pi*p.hbar^2*p.nSource/(p.gValleySource*p.mS); % Source chemical potential
p.Tc = 10;                       % Critical temperature in K
p.T  = 0.10;                     % Actual temperature in K
p.Delta0 = 1.764*p.kB*p.Tc;         % Superconducting gap
p.Delta = p.Delta0*tanh(1.74*sqrt(p.Tc/p.T - 1)); % Temperature dependent superconducting gap
deltaOverride=getenv('LOUISE_DELTA_MEV');
if ~isempty(deltaOverride)
    requestedDeltaMeV=str2double(deltaOverride);
    assert(isfinite(requestedDeltaMeV) && requestedDeltaMeV>0, ...
        'LOUISE_DELTA_MEV must be a positive finite value in meV.');
    targetDelta=requestedDeltaMeV*1e-3*p.e;
    % Keep the finite-temperature BCS gap relation self-consistent while
    % setting the requested Delta*. This changes Tc along the gap sweep.
    p.Tc=fzero(@(Tc)1.764*p.kB*Tc.*tanh(1.74*sqrt(Tc/p.T-1))-targetDelta, ...
        [p.T*(1+1e-9),max(2*p.T,2*targetDelta/(1.764*p.kB))]);
    p.Delta0=1.764*p.kB*p.Tc;
    p.Delta=p.Delta0*tanh(1.74*sqrt(p.Tc/p.T-1));
end
p.W = 1.0e-6;                     % Device width [m], applied only AFTER integration

%% Remote-impurity-limited mobility: mu_mob = mu_ref*(n/n_ref)^alpha
p.muMobRef = 1e4*1e-4;       % 10,000 cm^2/(V s) -> m^2/(V s)
p.nRef = 1e12*1e4;           % 1e12 cm^-2 -> m^-2
p.alphaMob = 1.0;            % Mobility-density exponent; distinct from alphaGeo
peakMobilityOverride=getenv('LOUISE_PEAK_MOBILITY_CM2VS');
if ~isempty(peakMobilityOverride)
    requestedPeakMobility=str2double(peakMobilityOverride);
    assert(isfinite(requestedPeakMobility) && requestedPeakMobility>0, ...
        'LOUISE_PEAK_MOBILITY_CM2VS must be positive and finite.');
    sourceDensityRatio=sheetDensityTF(p.muS,p)/p.nRef;
    p.muMobRef=(requestedPeakMobility/1e4)/sourceDensityRatio^p.alphaMob;
end
p.lambdaRef = p.hbar*p.muMobRef*sqrt(2*pi*p.nRef)/p.e; % Mean-free path
muMobMax_cm2Vs = 1e4*p.muMobRef * ...
    (sheetDensityTF(p.muS,p)/p.nRef)^p.alphaMob; % Mobility in cm^2/Vs

% Save this sweep, its figures, and rolling checkpoints separately.
scriptFullPath = which('SNS_2D_mobility_5placements_Rn_IcRn_gm_source_density_Louise');
assert(~isempty(scriptFullPath),'Cannot locate this _Louise script on the MATLAB path.');
scriptDir = fileparts(scriptFullPath);
runStamp = datestr(now,'yyyymmdd_HHMMSS');
runLabel='';
if ~isempty(peakMobilityOverride)
    runLabel=['_PeakMob_' strrep(peakMobilityOverride,'.','p')];
    if ~isempty(deltaOverride)
        runLabel=[runLabel '_Dstar_' strrep(deltaOverride,'.','p') 'meV']; % Induced superconducting gap
    end
elseif ~isempty(deltaOverride)
    runLabel=['_Delta_' strrep(deltaOverride,'.','p') 'meV'];
end
if ~isempty(toxOverride)
    runLabel=[runLabel '_tox_' strrep(toxOverride,'.','p') 'nm'];
end
runRoot=getenv('LOUISE_RUN_ROOT');
if isempty(runRoot), runRoot=scriptDir; end
if ~exist(runRoot,'dir'), mkdir(runRoot); end
opt.outputDir = fullfile(runRoot,['run_' runStamp runLabel]);
opt.figureDir = fullfile(opt.outputDir,'figures');
if ~exist(opt.figureDir,'dir'), mkdir(opt.figureDir); end
opt.checkpointEvery = 5;
diary(fullfile(opt.outputDir,'run.log'));

settingsLegend = {
    sprintf('L_g=%.0f nm; t_{ox}=%.0f nm; \\Delta=%.3f meV', ...
        p.Lg_nm,p.tox_nm,1e3*p.Delta/p.e)
    sprintf('\\mu_S=%.1f meV; max mobility=%.0f cm^2/(V s)', ...
        p.muS/p.e*1e3,muMobMax_cm2Vs)
    };

%% Numerical settings for the Gauss-Legendre calculation of the 2D current 
% from the 1D current intedral 
muRatio = linspace(1e-3*p.e/p.Delta,p.muS/p.Delta,51); % From 1 meV to muS % Chemical potential/Delta*
% Include exact CPR diagnostic points so those panels always contain data.
muRatio = unique(sort([muRatio,2,4,8,12]));
phi = linspace(0,pi,201); % Superconducting phase

opt.Npos = 5;                    % Number of random scattering positions
opt.seed = 20261009;            % Date
opt.scatterRegion = 'full';     % 'full' = [0,L], 'plateau' = [Ltrans,L-Ltrans]
opt.NxMob = 1601;                % Full-length integral for 1/lambda(x)
opt.Nslices = 80;               % BdG smooth-profile discretization
opt.NkyOpen = 24;               % GL (Gauss-Legendre) nodes for propagating k_y
opt.NkyEvan = 16;               % GL nodes for evanescent k_y
opt.NkyNormal = 96;             % Separate denser GL grid for G_N
opt.EvanescentCutoff = 8*p.Delta; % Evanescent mode cutoff
opt.omegaMax = 30*p.Delta; % Maximum matsubara frequency included
opt.minOmegaStop = 5*p.Delta;   % the code won't consider stopping the sum before reaching this energy.
opt.matsRelTol = 1e-6;          % the Matsubara sum must converge to a relative tolerance of 10^{-6}
opt.matsConsecutive = 6;        % the convergence condition must hold for six consecutive terms before the sum stops.
opt.Nmax = ceil(0.5*(opt.omegaMax/(pi*p.kB*p.T)-1 % converts the maximum energy and temperature into the highest Matsubara index
fprintf('Run output folder: %s\n',opt.outputDir);
fprintf('Sweep has %d mu_N/Delta points; rolling checkpoint every %d points.\n', ...
    numel(muRatio),opt.checkpointEvery);

% Normal conductance is calculated for propagating N transverse modes.
% Normal conductance is calculated only for calculating RN
% It includes contact mismatch and coherent effects of the single scatterer.
% Evanescent normal-state modes are omitted because a constant-transmission
% phenomenological scatterer is not uniquely defined there.

% Run options. Set false for faster preliminary runs.
opt.runSmooth = true;
opt.runAbrupt = true;
opt.runBallistic = true;

% Assumsings needed for code to run
assert(p.T>0 && p.T<p.Tc,'Require 0 < T < Tc.');
assert(p.Lungated_nm<100,'The underlap transition model requires Lungated < 100 nm.');
assert(p.Lplateau_nm>=0, ...
    'The cosine-ramp model requires Lg >= 2*lambda; this geometry has overlapping transitions.');
assert(all(muRatio > 0),'This script assumes mu_N > 0.');
assert(p.mS>0 && p.nSource>0 && p.gValleySource>=1, ...
    'Source effective mass, density and valley degeneracy must be positive.');

%% Five fixed random locations (common random numbers across entire sweep)
rng(opt.seed,'twister');
switch lower(opt.scatterRegion)
    case 'full'
        xScatter = sort(p.L*rand(1,opt.Npos));
    case 'plateau'
        xScatter = sort(p.Ltrans+(p.L-2*p.Ltrans)*rand(1,opt.Npos));
    otherwise
        error('scatterRegion must be full or plateau');
end
fprintf('\nScatterer positions [nm]: %s\n',mat2str(xScatter*1e9,4));
fprintf('Delta = %.5f meV, lambda_ref = %.3f um\n', ...
    p.Delta/p.e*1e3,p.lambdaRef*1e6);
fprintf(['Geometry: Lg=%.3f nm, Lungated=%.3f nm/side, tox=%.3f nm, ' ...
    'lambda=%.3f nm, Ltrans=%.3f nm/side, Lch=%.3f nm, plateau=%.3f nm\n'], ...
    p.Lg_nm,p.Lungated_nm,p.tox_nm,p.lambda_nm,p.Ltrans_nm,p.Lch_nm,p.Lplateau_nm);
fprintf('Source: nS = %.4g cm^-2, mS/m0 = %.4g, muS = %.4f meV\n', ...
    p.nSource_cm2,p.mS/p.m0,p.muS/p.e*1e3);
fprintf('N_mu=%d, N_ky=%d+%d, N_positions=%d\n', ...
    numel(muRatio),opt.NkyOpen,opt.NkyEvan,opt.Npos);

%% Gauss-Legendre nodes and weights
[xOpen,wOpen] = gaussLegendre(opt.NkyOpen);
[xEvan,wEvan] = gaussLegendre(opt.NkyEvan);
[xNormal,wNormal] = gaussLegendre(opt.NkyNormal);

%% Allocate output (critical-current densities in A/m)
Nmu = numel(muRatio); % Chemical potential/Delta* 
IcSmooth_mean = nan(1,Nmu);      % Average of five individually maximized Ic using smooth transision
IcAbrupt_mean = nan(1,Nmu);      % Average of five individually maximized Ic using abrupt transision
IcSmooth_ofMeanCPR = nan(1,Nmu); % Not used directly but takes average first, then maximizes 
IcAbrupt_ofMeanCPR = nan(1,Nmu); % Not used directly but takes average first, then maximizes 
IcIdeal = nan(1,Nmu); % Ideal Andreev current
IcSmooth_each = nan(opt.Npos,Nmu); % Smooth current, for all 5 placements
IcAbrupt_each = nan(opt.Npos,Nmu); % Abrupt current, for all 5 placements
tauSmooth = nan(1,Nmu); % Smooth transmission
tauAbrupt = nan(1,Nmu); % Abrupt transmission
rhoSmooth = nan(1,Nmu); % Dimensionless measure of accumulated scattering for smooth profile
rhoAbrupt = nan(1,Nmu); % Dimensionless measure of accumulated scattering for abrupt profile
nMid = nan(1,Nmu); % Carrier sheet density evaluated at that centre of the channel for each gate setting
lambdaMid = nan(1,Nmu); % Mean-free path evaluated at that centre of the channel for each gate setting

%% Normal-state two-terminal conductance G_N (conductance PER WIDTH, S/m)
% GN is only calculated to obtain RN
% The same fixed five scatterer positions are used as for Ic.
Gsm_each_perW = nan(opt.Npos,Nmu); % normal-state conductance per unit width for each scatterer placement for smooth profile
Gab_each_perW = nan(opt.Npos,Nmu); % normal-state conductance per unit width for each scatterer placement for abrupt profile
Gsm_mean_perW = nan(1,Nmu); % Averange normal-state conductance per unit width for smooth profile
Gab_mean_perW = nan(1,Nmu); % Averange normal-state conductance per unit width for abrupt profile
Gideal_perW = nan(1,Nmu); % Ballistic normal-state conductance per unit width
Gsm_incoh_perW = nan(1,Nmu);   % Gideal * tau_tot (no interface reflection). incoherent estimate of the normal conductance per unit width for smooth profile
% "Incoherent" means it applies that transmission as a simple probability factor, without calculating phase-coherent interference from the individual scatterer positions.
Gab_incoh_perW = nan(1,Nmu); % incoherent estimate of the normal conductance per unit width for abrupt profile

% Gate-voltage axis from midpoint Thomas-Fermi sheet density.
Vg = nan(1,Nmu);             % Vg(muN), offset p.Vfb
Cq = nan(1,Nmu);             % Differential channel capacitance [F/m^2]

% Select a few operating points for CPR and profile diagnostics.
muRatioCPR = [2,4,8,12];
CPR_smooth = cell(1,numel(muRatioCPR)); % stores the smooth-profile current–phase relations at the selected points
CPR_abrupt = cell(1,numel(muRatioCPR)); % stores the abrupt-profile current–phase relations at the selected points
CPR_ideal  = cell(1,numel(muRatioCPR)); % stores the current–phase relations at the selected points for the ideal Andreev

%% Gate-parameter sweep
for im = 1:Nmu
    muN = muRatio(im)*p.Delta;

    % Entire length, INCLUDING both smooth transitions.
    [tauSmooth(im),rhoSmooth(im)] = distributedTransmission( ...
        muN,'smooth',p,opt.NxMob);
    [tauAbrupt(im),rhoAbrupt(im)] = distributedTransmission( ...
        muN,'abrupt',p,opt.NxMob);

    nMid(im) = sheetDensityTF(muN,p); % Mid-channel sheet carrier density
    lambdaMid(im) = localMeanFreePath(nMid(im),p); % Mid-channel mean-free path

    % Gate potential = chemical-potential shift + oxide voltage.
    % This is a LOCAL midpoint approximation, not a self-consistent Poisson
    % solution for the entire contact-induced x-dependent potential.
    Vg(im) = p.Vfb + muN/p.e + p.e*nMid(im)/p.Cox; % Gate voltage from flat-band voltage offset, chemical potential
    % voltage drop across oxide
    Cq(im) = p.e^2*p.mN/(pi*p.hbar^2) * ...
        fermiOccupationAtBandEdge(muN,p); % Quantum capacitance per unit area

    % Split at k_y where longitudinal mu_N - E_y changes sign. This
    % improves Gauss-Legendre (GL) integration compared with one interval across threshold.
    % ky is conserved at each interface. Transverse energies differ
    % if mS ~= mN, so impose the S-lead cutoff using mS.
    kyF = sqrt(2*p.mN*muN)/p.hbar; % Largest transverse wave number that has a propagating mode in the normal channel, Fermi cutoff
    kyS = sqrt(2*p.mS*p.muS)/p.hbar; % Largest transverse wave number that has a propagating mode in the superconducting leads, Fermi cutoff
    kyMaxN = sqrt(2*p.mN*(muN+opt.EvanescentCutoff))/p.hbar; % Extends the normal-channel range beyond kyF by an energy opt.EvanescentCutoff. 
    % This includes some modes that are evanescent in the normal channel but can still affect the calculation.
    kyMax = min(kyMaxN,sqrt(0.95)*kyS); % Caps the integration range at the smaller of that extended normal-channel limit and sqrt(0.95)*kyS. 
    % The cap avoids integrating right up to the lead's transverse-mode cutoff.
    kyOpen = min(kyF,kyMax); % Marks the boundary between modes that propagate in the normal channel and evanescent modes
    [ky,wky] = transverseQuadrature(kyOpen,kyMax, ...
        xOpen,wOpen,xEvan,wEvan % The quadrature points ky and weights wky over the propagating and evanescent ranges, so the code can numerically integrate their contributions.

    Ismooth = zeros(opt.Npos,numel(phi)); % The CPR calculated for the smooth junction profile
    Iabrupt = zeros(opt.Npos,numel(phi)); % The CPR calculated for the abrupt junction profile
    Iideal = zeros(1,numel(phi)); % The CPR calculated for the the Andreev curve
    Gsm_each = zeros(opt.Npos,1);  % Conductance G/W for each placement [S/m] for smooth profile for each individual scattering position
    Gab_each = zeros(opt.Npos,1); % Conductance G/W for each placement [S/m] for abrupt profile for each individual scattering position
    Gquantum = 2*p.e^2/p.h;       % Spin-degenerate Landauer quantum [S]
    Gideal_perW(im) = Gquantum*min(kyF,kyS)/pi; % Ideal ballistic normal-state conductance per unit width
    Gsm_incoh_perW(im) = Gideal_perW(im)*tauSmooth(im); % Incoherent estimate of normal-state conductance per unit width for smooth profile
    Gab_incoh_perW(im) = Gideal_perW(im)*tauAbrupt(im); % Incoherent estimate of normal-state conductance per unit width for abrupt profile

    for ik = 1:numel(ky) % Go through all transverse modes with k vector ky
        EyN = p.hbar^2*ky(ik)^2/(2*p.mN); % Transverse kinetic energy of the mode with transverse wave number ky in normal region
        EyS = p.hbar^2*ky(ik)^2/(2*p.mS); % Transverse kinetic energy of the mode with transverse wave number ky in superconducting lead region
        muNx = muN-EyN; % Longitudinal energy available to a mode in the normal channel 
        muSxLead = p.muS-EyS;     % S reservoir (mass mS) Longitudinal energy available to a mode in the lead
        muSxRamp = p.muS-EyN;     % N-side ramp (mass mN) Longitudinal energy used for the superconducting-lead potential along the transition ramp
        weight = wky(ik)/pi;     % I/W = (1/pi)*integral_0^infty dk_y I_1D numerical integration weight
        if muSxLead<=0, continue; end % skips the current transverse mode if it has no positive longitudinal energy in the superconducting lead

        if opt.runSmooth % Calculate smooth CPR
            I1 = exactModeCPRs(phi,muNx,muSxLead,muSxRamp,tauSmooth(im), ...
                xScatter,'smooth',p,opt);
            Ismooth = Ismooth + weight*I1;
        end

        if opt.runAbrupt % Calculate abrupt CPR
            I1 = exactModeCPRs(phi,muNx,muSxLead,muSxRamp,tauAbrupt(im), ...
                xScatter,'abrupt',p,opt);
            Iabrupt = Iabrupt + weight*I1;
        end

        if opt.runBallistic && muNx>0 % Calculate ballistic CPR
            Iideal = Iideal + weight*idealModeCPR(phi,muNx,p,opt);
        end

    end

    % Dedicated, denser transverse quadrature (way to approximate an integral by adding up a carefully chosen set of sample values) for normal-state G_N.
    % No Matsubara sum needed: these transmissions are at real E=E_F.
    [kyNorm,wNorm] = mapGL(xNormal,wNormal,0,min(kyF,kyS)); % transverse wave numbers and qaudratue weights used to calculate the normal-state conductance
    for ik=1:numel(kyNorm) % go thourgh all transverse k vectors
        EyN = p.hbar^2*kyNorm(ik)^2/(2*p.mN); % transverse kinetic energy in normal region
        EyS = p.hbar^2*kyNorm(ik)^2/(2*p.mS); % transverse kinetic energy in supercnduing leads
        muNx = muN-EyN; % longitudinal energy in the normal channel
        muSxLead = p.muS-EyS; % longitudinal energy in the lead
        muSxRamp = p.muS-EyN; % energy available for motion along the transition ramp for the current transverse mode
        if muNx<=0 || muSxLead<=0, continue; end % skips the current kyNorm mode if it cannot propagate along the channel or the superconducting lead
        if opt.runSmooth
            Ts = normalModeTransmissions(muNx,muSxLead,muSxRamp, ...
                tauSmooth(im),xScatter,'smooth',p,opt); % normal-state transmission probability for the current
            Gsm_each = Gsm_each + (Gquantum/pi)*wNorm(ik)*Ts;
        end
        if opt.runAbrupt
            Ta = normalModeTransmissions(muNx,muSxLead,muSxRamp, ...
                tauAbrupt(im),xScatter,'abrupt',p,opt);
            Gab_each = Gab_each + (Gquantum/pi)*wNorm(ik)*Ta;
        end
    end

    if opt.runSmooth
        Gsm_each_perW(:,im) = Gsm_each;
        Gsm_mean_perW(im) = mean(Gsm_each);
        assert(Gsm_mean_perW(im) <= Gideal_perW(im)*(1+1e-5), ...
            'Smooth normal conductance exceeds ideal quantum limit.');
        IcSmooth_each(:,im) = max(abs(Ismooth),[],2);
        IcSmooth_mean(im) = mean(IcSmooth_each(:,im));
        IcSmooth_ofMeanCPR(im) = max(abs(mean(Ismooth,1)));
    end
    if opt.runAbrupt
        Gab_each_perW(:,im) = Gab_each;
        Gab_mean_perW(im) = mean(Gab_each);
        assert(Gab_mean_perW(im) <= Gideal_perW(im)*(1+1e-5), ...
            'Abrupt normal conductance exceeds ideal quantum limit.');
        IcAbrupt_each(:,im) = max(abs(Iabrupt),[],2);
        IcAbrupt_mean(im) = mean(IcAbrupt_each(:,im));
        IcAbrupt_ofMeanCPR(im) = max(abs(mean(Iabrupt,1)));
    end
    if opt.runBallistic
        IcIdeal(im) = max(abs(Iideal));
    end

    [d,idx] = min(abs(muRatioCPR-muRatio(im)));
    if d<1e-10
        CPR_smooth{idx} = Ismooth;
        CPR_abrupt{idx} = Iabrupt;
        CPR_ideal{idx} = Iideal;
    end

    fprintf(['mu/Delta=%5.2f  T_smooth=%.4f  T_abrupt=%.4f', ...
        '  Ic_smooth/W=%8.4g A/m\n'], ...
        muRatio(im),tauSmooth(im),tauAbrupt(im),IcSmooth_mean(im));
    if mod(im,opt.checkpointEvery)==0 || im==Nmu
        saveSweepProgress(opt.outputDir,opt.figureDir,settingsLegend, ...
            im,Nmu,muRatio,IcSmooth_mean,IcAbrupt_mean,IcIdeal, ...
            IcSmooth_each,IcAbrupt_each,tauSmooth,tauAbrupt, ...
            rhoSmooth,rhoAbrupt,Gsm_mean_perW,Gab_mean_perW, ...
            Gideal_perW,Vg,Cq,CPR_smooth,CPR_abrupt,CPR_ideal);
    end
end

%% Normal-state resistance and Ic*Rn
% Average conductance before inversion: R_N = 1/<G_N>.
% This is a width-dependent two-terminal resistance of the modeled junction,
% without additional external wiring resistance.
Rsm_ohm = 1./(p.W*Gsm_mean_perW);
Rab_ohm = 1./(p.W*Gab_mean_perW);
Rideal_ohm = 1./(p.W*Gideal_perW);
Rsm_incoh_ohm = 1./(p.W*Gsm_incoh_perW);
Rab_incoh_ohm = 1./(p.W*Gab_incoh_perW);

IcRn_sm_V = IcSmooth_mean./Gsm_mean_perW;
IcRn_ab_V = IcAbrupt_mean./Gab_mean_perW;
IcRn_ideal_V = IcIdeal./Gideal_perW;

% Alternative: average (Ic_j*Rn_j) realization by realization.
IcRn_sm_each_V = IcSmooth_each./Gsm_each_perW;
IcRn_ab_each_V = IcAbrupt_each./Gab_each_perW;
IcRn_sm_meanProduct_V = mean(IcRn_sm_each_V,1);
IcRn_ab_meanProduct_V = mean(IcRn_ab_each_V,1);

%% Josephson transconductance dIc/dVg [A/V]
% Vg uses n_mid(muN), and Ic already includes tau(muN).
% Numerically differentiated ensemble-averaged critical current.
dIc_dVg_sm_ApV = gradient(p.W*IcSmooth_mean,Vg);
dIc_dVg_ab_ApV = gradient(p.W*IcAbrupt_mean,Vg);
dIc_dVg_ideal_ApV = gradient(p.W*IcIdeal,Vg);
% The finite-grid numerical derivative may be noisy near FP resonances.
% Optional smoothing for presentation only (disabled by default).
% dIc_dVg_sm_plot = movmean(dIc_dVg_sm_ApV,3);

fprintf('Cox = %.3g F/m^2 (%.2f uF/cm^2); Cq = %.3g F/m^2\n', ...
    p.Cox,100*p.Cox,mean(Cq));
fprintf('d(muN/e)/dVg approx %.4f (midpoint gate lever arm)\n', ...
    mean(p.Cox./(p.Cox+Cq)));
fprintf('Ideal short-junction IcRn reference = %.4f mV\n', ...
    1e3*pi*p.Delta/p.e);

%% Plots
% Plotting, plotting, plotting...
% Mobility is largest in the superconducting source regions, where the
% local density reaches the density set by p.nSource_cm2.
muMobMax_cm2Vs = 1e4*p.muMobRef * ...
    (sheetDensityTF(p.muS,p)/p.nRef)^p.alphaMob;

figure('Color','w','Name','Critical current density'); hold on; grid on; box on;
if opt.runSmooth
    plot(muRatio,IcSmooth_mean,'LineWidth',2,'DisplayName', ...
        'Smooth: mean of 5 I_c');
end
if opt.runAbrupt
    plot(muRatio,IcAbrupt_mean,'LineWidth',1.8,'DisplayName', ...
        'Abrupt: mean of 5 I_c');
end
if opt.runBallistic
    plot(muRatio,IcIdeal,'k--','LineWidth',1.8, ...
        'DisplayName','Ideal ballistic Andreev');
end
xlabel('\mu_N/\Delta'); ylabel('I_c/W [A/m]');
title('2D SNS: position-averaged effective disorder');
lgd = legend('Location','best');
lgd.Title.String = {
    sprintf('L_g=%.0f nm; L_{ungated}=%.0f nm/side; L_{ch}=%.0f nm', ...
        p.Lg_nm,p.Lungated_nm,p.Lch_nm)
    sprintf('t_{ox}=%.0f nm; L_{trans}=%.1f nm/side; L_{plateau}=%.1f nm', ...
        p.tox_nm,p.Ltrans_nm,p.Lplateau_nm)
    sprintf('C_{ox}=%.3g uF/cm^2; \\Delta=%.3f meV', ...
        100*p.Cox,1e3*p.Delta/p.e)
    sprintf('\\mu_{mob,max}=%.3g cm^2/(V s)',muMobMax_cm2Vs)
    };

figure('Color','w','Name','Transmission and density');
tiledlayout(3,1);
nexttile; hold on; grid on;
plot(muRatio,tauSmooth,'LineWidth',2,'DisplayName','Smooth');
plot(muRatio,tauAbrupt,'--','LineWidth',1.6,'DisplayName','Abrupt');
ylabel('\tau_{tot}'); ylim([0 1]); legend('Location','best');
nexttile; grid on;
plot(muRatio,nMid/1e4,'LineWidth',1.8);
ylabel('n_{mid} [cm^{-2}]');
nexttile; grid on;
plot(muRatio,lambdaMid*1e9,'LineWidth',1.8);
ylabel('\lambda_{mid} [nm]'); xlabel('\mu_N/\Delta');

%% Normal conductance, on-resistance, IcRn and transconductance
figure('Color','w','Name','Normal-state Landauer conductance');
tiledlayout(2,1);
nexttile; hold on; grid on; box on;
plot(muRatio,p.W*Gsm_mean_perW*1e3,'LineWidth',2, ...
    'DisplayName','Smooth: full normal scattering');
plot(muRatio,p.W*Gab_mean_perW*1e3,'LineWidth',1.7, ...
    'DisplayName','Abrupt: full normal scattering');
plot(muRatio,p.W*Gideal_perW*1e3,'k--','LineWidth',1.7, ...
    'DisplayName','Ideal ballistic');
plot(muRatio,p.W*Gsm_incoh_perW*1e3,':','LineWidth',1.7, ...
    'DisplayName','Smooth: G_{ideal}\tau_{tot}');
ylabel('G_N [mS]'); title(sprintf('Normal conductance, W = %.1f \\mum',p.W*1e6));
legend('Location','best');
nexttile; hold on; grid on; box on;
plot(muRatio,Gsm_mean_perW./Gideal_perW,'LineWidth',2, ...
    'DisplayName','Smooth G_N/G_{ideal}');
plot(muRatio,Gab_mean_perW./Gideal_perW,'LineWidth',1.7, ...
    'DisplayName','Abrupt G_N/G_{ideal}');
plot(muRatio,tauSmooth,':','LineWidth',1.7, ...
    'DisplayName','Smooth \tau_{tot}');
ylabel('G_N/G_{ideal}'); xlabel('\mu_N/\Delta');
legend('Location','best');

figure('Color','w','Name','Normal-state on resistance'); hold on; grid on; box on;
plot(muRatio,Rsm_ohm,'LineWidth',2,'DisplayName','Smooth R_N');
plot(muRatio,Rab_ohm,'LineWidth',1.7,'DisplayName','Abrupt R_N');
plot(muRatio,Rideal_ohm,'k--','LineWidth',1.7, ...
    'DisplayName','Ideal ballistic R_N');
plot(muRatio,Rsm_incoh_ohm,':','LineWidth',1.7, ...
    'DisplayName','Smooth: [G_{ideal}\tau_{tot}]^{-1}');
xlabel('\mu_N/\Delta'); ylabel('R_N [\Omega]');
title(sprintf('Two-terminal on resistance, W=%.1f \\mum',p.W*1e6));
legend('Location','best');

figure('Color','w','Name','IcRn characteristic'); hold on; grid on; box on;
plot(muRatio,1e3*IcRn_sm_V,'LineWidth',2, ...
    'DisplayName','Smooth: I_c / mean(G_N)');
plot(muRatio,1e3*IcRn_ab_V,'LineWidth',1.7, ...
    'DisplayName','Abrupt: I_c / mean(G_N)');
plot(muRatio,1e3*IcRn_ideal_V,'k--','LineWidth',1.7, ...
    'DisplayName','Ideal finite-length Andreev');
yline(1e3*pi*p.Delta/p.e,'k:', ...
    'Ideal short-junction \pi\Delta/e','HandleVisibility','off');
xlabel('\mu_N/\Delta'); ylabel('I_cR_N [mV]');
title('Josephson I_cR_N product; R_N = 1/mean_j G_j');
legend('Location','best');

figure('Color','w','Name','Gate voltage and Josephson transconductance');
tiledlayout(2,1);
nexttile; grid on; box on;
plot(muRatio,Vg-p.Vfb,'LineWidth',1.8);
ylabel('V_g - V_{fb} [V]');
title(sprintf('Local gate mapping: C_{ox}=%.1f \\muF/cm^2',100*p.Cox));
nexttile; hold on; grid on; box on;
plot(Vg-p.Vfb,1e6*dIc_dVg_sm_ApV,'LineWidth',2, ...
    'DisplayName','Smooth');
plot(Vg-p.Vfb,1e6*dIc_dVg_ab_ApV,'LineWidth',1.7, ...
    'DisplayName','Abrupt');
plot(Vg-p.Vfb,1e6*dIc_dVg_ideal_ApV,'k--','LineWidth',1.7, ...
    'DisplayName','Ideal ballistic');
xlabel('V_g - V_{fb} [V]'); ylabel('dI_c/dV_g [\muA/V]');
legend('Location','best');

figure('Color','w','Name','Ic as a function of Vg'); hold on; grid on; box on;
plot(Vg-p.Vfb,1e6*p.W*IcSmooth_mean,'LineWidth',2,'DisplayName','Smooth');
plot(Vg-p.Vfb,1e6*p.W*IcAbrupt_mean,'LineWidth',1.7,'DisplayName','Abrupt');
plot(Vg-p.Vfb,1e6*p.W*IcIdeal,'k--','LineWidth',1.7, ...
    'DisplayName','Ideal ballistic');
xlabel('V_g - V_{fb} [V]'); ylabel('I_c [\muA]');
legend('Location','best');

figure('Color','w','Name','Five realizations'); hold on; grid on;
for j=1:opt.Npos
    plot(muRatio,IcSmooth_each(j,:),'LineWidth',1.0, ...
        'DisplayName',sprintf('x_s=%.1f nm',xScatter(j)*1e9));
end
plot(muRatio,IcSmooth_mean,'k-','LineWidth',2.2, ...
    'DisplayName','Mean I_c');
xlabel('\mu_N/\Delta'); ylabel('I_c/W [A/m]');
title('Smooth profile: five fixed random positions'); legend('Location','best');

figure('Color','w','Name','Average vs maximize'); hold on; grid on;
plot(muRatio,IcSmooth_mean,'LineWidth',2, ...
    'DisplayName','Mean of 5 critical currents');
plot(muRatio,IcSmooth_ofMeanCPR,'--','LineWidth',1.8, ...
    'DisplayName','Critical current of average CPR');
xlabel('\mu_N/\Delta'); ylabel('I_c/W [A/m]');
legend('Location','best');

figure('Color','w','Name','Selected current-phase relations');
tiledlayout(2,2);
for j=1:numel(muRatioCPR)
    nexttile; hold on; grid on;
    if ~isempty(CPR_smooth{j})
        plot(phi/pi,mean(CPR_smooth{j},1),'LineWidth',2, ...
            'DisplayName','Smooth, mean CPR');
    end
    if ~isempty(CPR_abrupt{j})
        plot(phi/pi,mean(CPR_abrupt{j},1),'LineWidth',1.5, ...
            'DisplayName','Abrupt, mean CPR');
    end
    if ~isempty(CPR_ideal{j})
        plot(phi/pi,CPR_ideal{j},'k--','LineWidth',1.4, ...
            'DisplayName','Ideal Andreev');
    end
    xlabel('\phi/\pi'); ylabel('I/W [A/m]');
    title(sprintf('\\mu_N/\\Delta=%g',muRatioCPR(j)));
    if j==1, legend('Location','best'); end
end

figure('Color','w','Name','Local mobility profile');
xplot = linspace(0,p.L,701);
exampleRatio = 4;
muX = channelChemicalPotential(xplot,exampleRatio*p.Delta,p.muS,'smooth',p);
nX = sheetDensityTF(muX,p);
muMobX = p.muMobRef*(nX/p.nRef).^p.alphaMob;
lamX = localMeanFreePath(nX,p);
tiledlayout(3,1);
nexttile; plot(xplot*1e9,muX/p.Delta,'LineWidth',1.8); grid on;
ylabel('\mu(x)/\Delta'); title(sprintf('Smooth profile, \\mu_N/\\Delta=%g',exampleRatio));
nexttile; plot(xplot*1e9,muMobX*1e4,'LineWidth',1.8); grid on;
ylabel('\mu_{mob}(x) [cm^2/Vs]');
nexttile; plot(xplot*1e9,lamX*1e9,'LineWidth',1.8); grid on;
ylabel('\lambda(x) [nm]'); xlabel('x [nm]');

% Add shared settings to every plotted axes legend, then export all figures.
figs=findall(groot,'Type','figure');
for jf=1:numel(figs)
    axs=findall(figs(jf),'Type','axes');
    for ja=1:numel(axs)
        if isempty(findall(axs(ja),'Type','line')), continue; end
        lg=findall(figs(jf),'Type','legend');
        lg=lg(arrayfun(@(x) isequal(x.Axes,axs(ja)),lg));
        if isempty(lg)
            h=findall(axs(ja),'Type','line','-not','HandleVisibility','off');
            if isempty(h), continue; end
            lg=legend(axs(ja),h(1),'Location','best');
        end
        lg(1).Title.String=settingsLegend;
    end
end
exportRunFigures(opt.figureDir);
%% Save the computed quantities for subsequent Vg modeling
results = struct();
results.muRatio = muRatio;
results.muN_J = muRatio*p.Delta;
results.xScatter_m = xScatter;
results.phi = phi;
results.tauSmooth = tauSmooth;
results.tauAbrupt = tauAbrupt;
results.rhoSmooth = rhoSmooth;
results.rhoAbrupt = rhoAbrupt;
results.nMid_m2 = nMid;
results.lambdaMid_m = lambdaMid;
results.Cox_Fpm2 = p.Cox;
results.geometry = struct('Lg_nm',p.Lg_nm,'Lungated_nm',p.Lungated_nm, ...
    'tox_nm',p.tox_nm,'tch_nm',p.tch_nm,'alphaGeo',p.alphaGeo, ...
    'epsilonChannelRelative',p.epsilonChannelRelative, ...
    'epsilonOxRelative',p.epsilonOxRelative,'lambda_nm',p.lambda_nm, ...
    'Ltrans_nm',p.Ltrans_nm,'Lch_nm',p.Lch_nm, ...
    'Lplateau_nm',p.Lplateau_nm);
results.mobilityModel = struct('muMobRef_cm2Vs',1e4*p.muMobRef, ...
    'nRef_cm2',p.nRef/1e4,'alphaMob',p.alphaMob, ...
    'alphaGeo',p.alphaGeo,'muMobMax_cm2Vs',muMobMax_cm2Vs);
results.nSource_m2 = p.nSource;
results.muSource_J = p.muS;
results.mSource_kg = p.mS;
results.Vg_V = Vg;
results.Cq_Fpm2 = Cq;
results.Gsmooth_each_perW_Spm = Gsm_each_perW;
results.Gabrupt_each_perW_Spm = Gab_each_perW;
results.Gsmooth_mean_perW_Spm = Gsm_mean_perW;
results.Gabrupt_mean_perW_Spm = Gab_mean_perW;
results.Gideal_perW_Spm = Gideal_perW;
results.Gsmooth_incoh_perW_Spm = Gsm_incoh_perW;
results.Gabrupt_incoh_perW_Spm = Gab_incoh_perW;
results.Rsmooth_ohm = Rsm_ohm;
results.Rabrupt_ohm = Rab_ohm;
results.Rideal_ohm = Rideal_ohm;
results.Rsmooth_incoh_ohm = Rsm_incoh_ohm;
results.Rabrupt_incoh_ohm = Rab_incoh_ohm;
results.IcRn_smooth_V = IcRn_sm_V;
results.IcRn_abrupt_V = IcRn_ab_V;
results.IcRn_ideal_V = IcRn_ideal_V;
results.IcRn_smooth_meanProduct_V = IcRn_sm_meanProduct_V;
results.IcRn_abrupt_meanProduct_V = IcRn_ab_meanProduct_V;
results.dIc_dVg_smooth_ApV = dIc_dVg_sm_ApV;
results.dIc_dVg_abrupt_ApV = dIc_dVg_ab_ApV;
results.dIc_dVg_ideal_ApV = dIc_dVg_ideal_ApV;
results.IcSmooth_mean_Apm = IcSmooth_mean;
results.IcAbrupt_mean_Apm = IcAbrupt_mean;
results.IcSmooth_each_Apm = IcSmooth_each;
results.IcAbrupt_each_Apm = IcAbrupt_each;
results.IcSmooth_ofMeanCPR_Apm = IcSmooth_ofMeanCPR;
results.IcAbrupt_ofMeanCPR_Apm = IcAbrupt_ofMeanCPR;
results.IcIdeal_Apm = IcIdeal;
results.IcSmooth_total_A = p.W*IcSmooth_mean;
results.IcAbrupt_total_A = p.W*IcAbrupt_mean;
results.p = p;
results.opt = opt;
results.outputDir = opt.outputDir;
save(fullfile(opt.outputDir,'SNS_2D_Rn_IcRn_gm_results.mat'),'results','-v7.3');

%% Optional easy-to-inspect numerical table (also saved as CSV)
summaryTable = table(muRatio(:), (Vg(:)-p.Vfb), ...
    tauSmooth(:), tauAbrupt(:), ...
    (p.W*Gsm_mean_perW(:)*1e3), (p.W*Gab_mean_perW(:)*1e3), ...
    (p.W*Gideal_perW(:)*1e3), ...
    Rsm_ohm(:), Rab_ohm(:), Rideal_ohm(:), ...
    (p.W*IcSmooth_mean(:)*1e6), (p.W*IcAbrupt_mean(:)*1e6), ...
    (p.W*IcIdeal(:)*1e6), ...
    (1e3*IcRn_sm_V(:)), (1e3*IcRn_ab_V(:)), (1e3*IcRn_ideal_V(:)), ...
    (1e6*dIc_dVg_sm_ApV(:)), (1e6*dIc_dVg_ab_ApV(:)), ...
    'VariableNames',{'muOverDelta','VgRelative_V', ...
    'tauSmooth','tauAbrupt','Gsmooth_mS','Gabrupt_mS','Gideal_mS', ...
    'Rsmooth_Ohm','Rabrupt_Ohm','Rideal_Ohm', ...
    'IcSmooth_uA','IcAbrupt_uA','IcIdeal_uA', ...
    'IcRnSmooth_mV','IcRnAbrupt_mV','IcRnIdeal_mV', ...
    'gmSmooth_uApV','gmAbrupt_uApV'});
writetable(summaryTable,fullfile(opt.outputDir,'SNS_2D_Rn_IcRn_gm_summary.csv'));
disp(summaryTable(1:5:end,:));

%% ===================== LOCAL FUNCTIONS ===============================

function saveSweepProgress(outDir,figDir,settingsLegend,im,Nmu,muRatio,Ics,Ica,Ici, ...
    IcsEach,IcaEach,tauS,tauA,rhoS,rhoA,Gs,Ga,Gi,Vg,Cq,CPRs,CPRa,CPRi)
progress=struct('completedIndex',im,'totalPoints',Nmu,'muRatio',muRatio(1:im), ...
    'IcSmooth_mean_Apm',Ics(1:im),'IcAbrupt_mean_Apm',Ica(1:im), ...
    'IcIdeal_Apm',Ici(1:im),'IcSmooth_each_Apm',IcsEach(:,1:im), ...
    'IcAbrupt_each_Apm',IcaEach(:,1:im),'tauSmooth',tauS(1:im), ...
    'tauAbrupt',tauA(1:im),'rhoSmooth',rhoS(1:im),'rhoAbrupt',rhoA(1:im), ...
    'Gsmooth_mean_perW_Spm',Gs(1:im),'Gabrupt_mean_perW_Spm',Ga(1:im), ...
    'Gideal_perW_Spm',Gi(1:im),'Vg_V',Vg(1:im),'Cq_Fpm2',Cq(1:im), ...
    'CPR_smooth',{CPRs},'CPR_abrupt',{CPRa},'CPR_ideal',{CPRi}, ...
    'settings',{settingsLegend});
tmp=fullfile(outDir,'sweep_checkpoint_tmp.mat');
save(tmp,'progress','-v7.3'); movefile(tmp,fullfile(outDir,'sweep_checkpoint.mat'),'f');
% Headless batch sweeps can encounter MATLAB graphics handshaking timeouts.
% The numerical checkpoint is already safely written, so omit this optional
% progress image when the runner requests a graphics-free calculation.
if strcmp(getenv('LOUISE_SKIP_RUN_FIGURES'),'1'), return; end
f=figure('Visible','off','Color','w','Name','Sweep progress');
tiledlayout(f,2,1);
ax=nexttile; hold(ax,'on'); grid(ax,'on');
plot(ax,muRatio(1:im),Ics(1:im),'DisplayName','Smooth');
plot(ax,muRatio(1:im),Ica(1:im),'DisplayName','Abrupt');
plot(ax,muRatio(1:im),Ici(1:im),'k--','DisplayName','Ideal');
xlabel(ax,'\mu_N/\Delta'); ylabel(ax,'I_c/W [A/m]');
title(ax,sprintf('Completed %d of %d points',im,Nmu)); lg=legend(ax,'Location','best'); lg.Title.String=settingsLegend;
ax=nexttile; hold(ax,'on'); grid(ax,'on');
plot(ax,muRatio(1:im),tauS(1:im),'DisplayName','Smooth');
plot(ax,muRatio(1:im),tauA(1:im),'--','DisplayName','Abrupt');
xlabel(ax,'\mu_N/\Delta'); ylabel(ax,'\tau_{tot}'); ylim(ax,[0 1]);
exportgraphics(f,fullfile(figDir,'sweep_progress.png'),'Resolution',150); close(f);
end

function exportRunFigures(figDir)
% Batch sweep results are consumed by the separate manuscript figure builder.
% Avoid exporting many diagnostic windows inside the long numerical job.
if strcmp(getenv('LOUISE_SKIP_RUN_FIGURES'),'1')
    fprintf('Skipping per-run figure export (LOUISE_SKIP_RUN_FIGURES=1).\n');
    return;
end

figs=findall(groot,'Type','figure');
for k=1:numel(figs)
    nm=get(figs(k),'Name'); if isempty(nm), nm=sprintf('figure_%02d',k); end
    nm=regexprep(nm,'[^A-Za-z0-9_-]','_');
    exportgraphics(figs(k),fullfile(figDir,[nm '.png']),'Resolution',200);
    exportgraphics(figs(k),fullfile(figDir,[nm '.pdf']),'ContentType','vector');
end
end

% Function which calculates sheet carrier density
function n = sheetDensityTF(mu,p)
% Single 2D subband, spin degeneracy 2, with stable finite-T log(1+exp(x)).
    u = mu/(p.kB*p.T); % Expresses the chemical potential in units of thermal energy \(k_BT\)
    logF = max(u,0) + log1p(exp(-abs(u))); % evaluates ln(1+e^u), but in a numerically stable way. Directly calculating e^u can overflow when u is large
    n = (p.mN*p.kB*p.T/(pi*p.hbar^2))*logF;  % Sheet carrier density in m^-2
end

% Function which calcualtes the mean-free path
function lambda = localMeanFreePath(n,p)
    lambda = p.lambdaRef*(n/p.nRef).^(p.alphaMob+0.5); % Mean-free path
end

% Function which calculates effective transmission
function [tau,rho] = distributedTransmission(muN,profile,p,Nx)
% Incoherent series resistance across 0 <= x <= L (including ramps).
    x = linspace(0,p.L,Nx); % Position in x-direction
    muX = channelChemicalPotential(x,muN,p.muS,profile,p); % Local chemical potential as a function of position x
    nx = sheetDensityTF(muX,p); % Local sheet carrier density as a function of position x
    lambda = localMeanFreePath(nx,p); % Local mean-free path as a function of position x
    rho = trapz(x,1./lambda); % Dimensionless measure of accumulated scattering: a larger value means more scattering
    tau = 1/(1+rho); % Converts that accumulated scattering into an effective transmission between 0 and 1
end

% Function that calculates the local chemical potential mu(x) along the junction
function mu = channelChemicalPotential(x,muN,muS,profile,p)
% x may be scalar or array; mu is relative to the Fermi level.
    if strcmpi(profile,'abrupt')
        % Abrupt reference: muS in both ungated underlaps, muN under the gate.
        mu = muS+zeros(size(x));
        gate = x>=p.Lungated_nm*1e-9 & ...
            x<=(p.Lungated_nm+p.Lg_nm)*1e-9;
        mu(gate) = muN;
        return
    end
    if ~strcmpi(profile,'smooth'), error('Bad profile: %s',profile); end
    mu = muN+zeros(size(x));
    left = (x<p.Ltrans);
    right = (x>p.L-p.Ltrans);
    mu(left) = muN+0.5*(muS-muN).* ...
        (1+cos(pi*x(left)/p.Ltrans));
    xr = p.L-x(right);
    mu(right) = muN+0.5*(muS-muN).* ...
        (1+cos(pi*xr/p.Ltrans));
end

% Function that creates the transverse-mode sampling points and their
% integration weights for the Gauss-Legendre calculation
function [ky,wky] = transverseQuadrature(kyF,kyMax,xO,wO,xE,wE)
% Integrate separately over propagating and evanescent transverse modes.
% x0 are the standard Gauss–Legendre nodes on [-1,1]
% w0 are their integration weights
    [ky1,w1] = mapGL(xO,wO,0,min(kyF,kyMax));
    [ky2,w2] = mapGL(xE,wE,min(kyF,kyMax),kyMax);
    ky = [ky1;ky2];
    wky = [w1;w2];
end

% Fucntion that takes Gauss–Legendre sample points and weights defined on the standard interval [-1,1], 
% then maps them onto the interval [a,b]
function [kn,wn] = mapGL(x,w,a,b)
    if b<=a
        kn = zeros(0,1); wn = zeros(0,1); return;
    end
    kn = a+0.5*(b-a)*(x+1);
    wn = 0.5*(b-a)*w;
end

% Function that calculates the current–phase relation for one transverse mode
function Ipos = exactModeCPRs(phi,muN,muSlead,muSramp,tau,positions,profile,p,opt)
% Single transverse mode: return Npos-by-Nphi CPRs [A], one per position.
% The electron and hole normal-region transfer matrices are conjugates
% because the normal electrostatic potential / scatterer is nonmagnetic.
    NP = numel(positions); % number of scatterer positions
    Ipos = zeros(NP,numel(phi)); % allocate memory
    nSmall = 0; % counter
    cosphi = cos(phi); sinphi = sin(phi);

    % creates the ordered positions that divide the junction into segments
    % for the transfer-matrix calculation for both the smooth and abrupt
    % case
    if strcmpi(profile,'abrupt')
        xEdge = [0,sort(positions),p.L];
    else
        xEdge = sort([linspace(0,p.L,opt.Nslices+1),positions]);
    end

    % Each random position occurs as an exact slice boundary.
    edgeIndex = zeros(1,NP);
    for j=1:NP
        [~,edgeIndex(j)] = min(abs(xEdge-positions(j)));
    end
    dx = diff(xEdge);
    muSeg = channelChemicalPotential( ...
        0.5*(xEdge(1:end-1)+xEdge(2:end)),muN,muSramp,profile,p); % local chemical-potential one for each segment
    muAtScatter = channelChemicalPotential(positions,muN,muSramp,profile,p); % local chemical-potential at each scatterer position
    Nseg = numel(dx); % number of spacial segments

    for n=0:opt.Nmax % goes through all matsubara frequcies
        omega = (2*n+1)*pi*p.kB*p.T; % matsubara energy
        prefix = zeros(2,2,Nseg+1); % stores the cumulative electron-propagation transfer matrices from the left end of the channel to each segment boundary
        suffix = zeros(2,2,Nseg+1); % stores the cumulative electron-propagation transfer matrices from each segment boundary to the right end of the channel
        prefix(:,:,1) = eye(2); % no propagation has been accumulated yet
        suffix(:,:,Nseg+1) = eye(2); % at that boundary there are no segments left to propagate through
        
        for j=1:Nseg % goes fhourgh all spacial segments
            Pslice = electronPropagation(muSeg(j),omega,dx(j),p); % 2×2 electron-propagation transfer matrix for one channel segment
            prefix(:,:,j+1) = normalizeTransferMatrix( ...
                Pslice*prefix(:,:,j));
        end

        for j=Nseg:-1:1
            Pslice = electronPropagation(muSeg(j),omega,dx(j),p);
            suffix(:,:,j) = normalizeTransferMatrix( ...
                suffix(:,:,j+1)*Pslice);
        end

        dI = zeros(NP,numel(phi));
        for j=1:NP
            PeL = prefix(:,:,edgeIndex(j));
            PeR = suffix(:,:,edgeIndex(j));
            Be = pointScattererElectron(muAtScatter(j),omega,tau,p);
            % The secular CPR depends on determinant ratios, so any common
            % nonzero scalar multiplying Pe cancels. Normalize this product
            % to prevent overflow/underflow for long depleted channels.
            Pe = normalizeTransferMatrix(PeR*Be*PeL);
            [A,B] = phaseCoefficients(Pe,omega,muSlead,p);
            % Derivative of ln[A+B*cos(phi)] with no finite difference.
            scale = max(abs([A,B]));
            if scale==0 || ~isfinite(scale)
                error('Invalid secular determinant at muN=%g, omega=%g.',muN,omega);
            end
            A=A/scale; B=B/scale;
            dI(j,:) = (4*p.e*p.kB*p.T/p.hbar)* ...
                real(B*sinphi./(A+B*cosphi));
        end
        Ipos = Ipos+dI;

        if omega>opt.minOmegaStop
            magnitude = max(abs(dI),[],2);
            scale = max(max(abs(Ipos),[],2),1e-6*p.e*p.Delta/p.hbar);
            if all(magnitude<opt.matsRelTol*scale)
                nSmall=nSmall+1;
            else
                nSmall=0;
            end
            if nSmall>=opt.matsConsecutive, break; end
        end
    end
end

function M = normalizeTransferMatrix(M)
% Remove a common scalar from a transfer matrix. This preserves the CPR,
% whose phase derivative uses only ratios of secular determinants.
    s = max(abs(M(:)));
    if ~isfinite(s) || s==0
        error('Transfer matrix became zero or nonfinite during propagation.');
    end
    M = M/s;
end

% Function that calculates the 2×2 transfer matrix for an electron moving through one uniform channel segment
function Pe = electronPropagation(mu,omega,dx,p)
% State is [psi; L*psi']. Dimensionless and well-conditioned in SI units.
    k = sqrt(2*p.mN*(mu+1i*omega))/p.hbar; % complex electron wave number for that segment at the Matsubara energy omega
    z = k*p.L; % dimensionless wave number
    theta = k*dx; % complex phase accumulated across one segment
    c = cos(theta);
    % calculates the off-diagonal element used in the propagation matrix while avoiding division by a nearly zero number
    if abs(z)<1e-12
        off = dx/p.L;
    else
        off = sin(theta)/z;
    end
    Pe = [c,off; -z*sin(theta),c];
end

function Be = pointScattererElectron(mu,omega,tau,p)
% Fixed normal-state t=sqrt(tau), r=i*sqrt(1-tau).
% Exact amplitude-to-state similarity transformation (NO matrix inverse).
% In [psi; L psi'] the point interface is a real SL(2,R) matrix at E real.
    if tau>=1-1e-14
        Be=eye(2); return;
    end
    if tau<=0, error('The scatterer needs 0 < tau <= 1'); end
    z = p.L*sqrt(2*p.mN*(mu+1i*omega))/p.hbar;
    s = sqrt(1-tau);
    Be = [1,-s/z; -s*z,1]/sqrt(tau);
    % Particle-hole counterpart is conj(Be) at E=i*omega, because the
    % boundary condition in the physical (psi,psi') basis is real at real E.
end

function [A,B] = phaseCoefficients(Pe,omega,muS,p)
% Exact 4x4 BdG determinant evaluated at 0 and pi; D=A+B*cos(phi).
    PN = zeros(4,4);
    PN([1,3],[1,3]) = Pe;
    PN([2,4],[2,4]) = conj(Pe);
    Omega = sqrt(omega^2+p.Delta^2);
    kp = sqrt(2*p.mS*(muS+1i*Omega))/p.hbar;
    km = sqrt(2*p.mS*(muS-1i*Omega))/p.hbar;
    if imag(kp)<0, kp=-kp; end
    if imag(km)<0, km=-km; end
    gp=1i*(omega-Omega)/p.Delta;
    gm=1i*(omega+Omega)/p.Delta;
    kpf = 1i*kp*p.L*p.mN/p.mS;
    kmf = 1i*km*p.L*p.mN/p.mS;
    D0  = leadDet(PN,0,gp,gm,kpf,kmf);
    Dpi = leadDet(PN,pi,gp,gm,kpf,kmf);
    A=0.5*(D0+Dpi);
    B=0.5*(D0-Dpi);
end

function D = leadDet(PN,phi,gp,gm,kpf,kmf)
    thetaL=-phi/2; thetaR=phi/2;
    gLp=gp*exp(-1i*thetaL); gLm=gm*exp(-1i*thetaL);
    gRp=gp*exp(-1i*thetaR); gRm=gm*exp(-1i*thetaR);
    BL=[1,1;gLp,gLm;-kpf,-kmf;-kpf*gLp,-kmf*gLm];
    BR=[1,1;gRp,gRm; kpf, kmf; kpf*gRp, kmf*gRm];
    M=[PN*BL,-BR];
    cn=sqrt(sum(abs(M).^2,1));
    M=bsxfun(@rdivide,M,cn);
    D=det(M);
end

% Ideal Andreev
function I = idealModeCPR(phi,muNx,p,opt)
% Perfect Andreev finite-length 1D reference, spin degeneracy included.
    if muNx<=0, I=zeros(size(phi)); return; end % If the mode has no positive longitudinal energy in the normal channel 
    % then assigns it zero current at every phase and skips the rest of the calculation
    vF=sqrt(2*muNx/p.mN); % Fermi velocity
    omega=(2*(0:opt.Nmax)+1)*pi*p.kB*p.T; % list of Matsubara energies
    a=2*asinh(omega/p.Delta)+2*omega*p.L/(p.hbar*vF); % The first term accounts for the superconducting gap, and the second accounts for propagation across a junction of length L at velocity v_F
    qq=exp(-a(:));
    c=cos(phi); s=sin(phi);
    I=(4*p.e*p.kB*p.T/p.hbar)*sum( ...
        bsxfun(@rdivide,2*qq*s, ...
        bsxfun(@plus,1+qq.^2,2*qq*c)),1);
end

% Fermi-dirac distribution
function f = fermiOccupationAtBandEdge(mu,p)
% Stable f(mu) = 1/(1+exp[-mu/(kB*T)]), used in Cq.
    t = mu/(p.kB*p.T);
    if t>=0
        f=1/(1+exp(-t));
    else
        y=exp(t); f=y/(1+y);
    end
end

function T = normalModeTransmissions(muNx,muSxLead,muSxRamp,tau,positions,profile,p,opt)
% Normal-state spin-degenerate Landauer transmission at E=E_F, one ky mode.
% Returns a transmission probability for each scatterer position.
%
% This evaluates the same real-potential scattering problem as the BdG
% normal region. It uses propagating lead and channel modes (muNx>0).
% For evanescent modes in the center, the phenomenological fixed-tau
% scatterer has no unique physical continuation; they are omitted here.
    assert(muNx>0 && muSxLead>0);
    NP = numel(positions);
    T = zeros(NP,1);

    if strcmpi(profile,'abrupt')
        xEdge = [0,sort(positions),p.L];
    else
        xEdge = sort([linspace(0,p.L,opt.Nslices+1),positions]);
    end
    edgeIndex = zeros(1,NP);
    for j=1:NP
        [~,edgeIndex(j)] = min(abs(xEdge-positions(j)));
    end
    dx = diff(xEdge);
    muSeg = channelChemicalPotential( ...
        0.5*(xEdge(1:end-1)+xEdge(2:end)),muNx,muSxRamp,profile,p);
    muAtScatter = channelChemicalPotential(positions,muNx,muSxRamp,profile,p);
    Nseg = numel(dx);

    prefix = zeros(2,2,Nseg+1);
    suffix = zeros(2,2,Nseg+1);
    prefix(:,:,1)=eye(2);
    suffix(:,:,Nseg+1)=eye(2);
    for i=1:Nseg
        U = electronPropagation(muSeg(i),0,dx(i),p);
        prefix(:,:,i+1) = U*prefix(:,:,i);
    end
    for i=Nseg:-1:1
        U = electronPropagation(muSeg(i),0,dx(i),p);
        suffix(:,:,i) = suffix(:,:,i+1)*U;
    end

    % State basis is [psi; L*(mN/m)*psi']. In the superconducting-normal
    % reservoirs at Delta=0, the normalized plane-wave derivative is
    % +/- i*zS (dimensionless), with zS including effective-mass ratio.
    kS = sqrt(2*p.mS*muSxLead)/p.hbar;
    zS = p.L*(p.mN/p.mS)*kS;
    for j=1:NP
        Pe = suffix(:,:,edgeIndex(j)) * ...
            pointScattererElectron(muAtScatter(j),0,tau,p) * ...
            prefix(:,:,edgeIndex(j));
        % C^(-1)*Pe*C, C=[1,1; i*zS,-i*zS].
        % M22=(a+d+i*(c/zS-zS*b))/2, so T=|1/M22|^2.
        denominator = 0.5*(Pe(1,1)+Pe(2,2) + ...
            1i*(Pe(2,1)/zS-zS*Pe(1,2)));
        T(j) = abs(1/denominator)^2;
        if ~isfinite(T(j)) || T(j)<0 || T(j)>1+1e-6
            error('Non-unitary normal transmission: muNx=%.5g, T=%g.', ...
                muNx,T(j));
        end
        T(j) = min(1,T(j));
    end
end

% Gauss-legendre is a numerical method to calculate the integral over 
% the current from 1D to 2D
function [x,w] = gaussLegendre(N) 
    n=1:N-1;
    beta=n./sqrt(4*n.^2-1);
    [V,D]=eig(diag(beta,1)+diag(beta,-1));
    [x,ind]=sort(diag(D));
    V=V(:,ind);
    w=2*(V(1,:).^2).';
end
