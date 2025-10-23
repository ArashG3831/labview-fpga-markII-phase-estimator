function final_v7()
% FPGA vs MATLAB with decimated logging (e.g., keep 1 of every DECIM samples)

%% ---------------- user settings ----------------
phi_file   = 'phi.txt';
sgphi_file = 'sgphi.txt';
sgi_file   = 'sgi.txt';
noise_file = 'noise.txt';

FS_EXPECT   = 125e6;      % nominal FPGA rate (pre-decimation)
FS_TOL_PPM  = 2e3;        % warn if off by >0.2%

N_full      = 6250;       % samples per pulse on FPGA (full-rate)
DECIM       = 5;          % you logged every 5th sample
mu          = 5.0;        % loop gain
A           = 1.0;        % homodyne slope
EPSG        = 1e-12;
TRIM_FRAC   = 0.01;

% per-pulse ψ_true sequence (repeats): 0.3→0.2→0.1→…
PSI_SEQ     = [0.3, 0.2, 0.1];

% Delay handling (in **full-rate** samples)
L_MODE   = 'fixed';       % 'fixed' or 'grid'
L_FIXED  = 23;            % known HW pipeline delay (full-rate)
L_GRID   = 0:80;          % grid (full-rate), if L_MODE='grid'

FIT_NOISE_SCALE = true;   % fit α (counts→units)

%% ---------------- load & align ----------------
[t_phi,   phi_hw,   Fs1] = load_trace(phi_file);
[t_sgphi, Sgphi_hw, Fs2] = load_trace(sgphi_file);
[t_sgi,   Sgi_hw,   Fs3] = load_trace(sgi_file);
[t_n,     n_counts, Fs4] = load_trace(noise_file);

Fs_dec = median([Fs1,Fs2,Fs3,Fs4]);          % this is the post-decimation Fs
ppm = 1e6*abs(Fs_dec - FS_EXPECT/DECIM)/(FS_EXPECT/DECIM);
if ppm > FS_TOL_PPM
    warning('Fs(off by %.0f ppm): Fs_dec=%.6f MHz expected=%.6f MHz.', ...
            ppm, Fs_dec/1e6, (FS_EXPECT/DECIM)/1e6);
end

% integer-sample start alignment on the decimated grid
t0 = min([t_phi(1), t_sgphi(1), t_sgi(1), t_n(1)]);
ishift = @(t) round((t(1)-t0)*Fs_dec);
phi_hw   = circshift(phi_hw,   -ishift(t_phi));
Sgphi_hw = circshift(Sgphi_hw, -ishift(t_sgphi));
Sgi_hw   = circshift(Sgi_hw,   -ishift(t_sgi));
n_counts = circshift(n_counts, -ishift(t_n));

% equal length + trim
M = min([numel(phi_hw), numel(Sgphi_hw), numel(Sgi_hw), numel(n_counts)]);
phi_hw   = phi_hw(1:M);  Sgphi_hw = Sgphi_hw(1:M);  Sgi_hw = Sgi_hw(1:M);  n_counts = n_counts(1:M);
t  = (0:M-1)'/Fs_dec;
edge = max(64, round(TRIM_FRAC*M));
keep = (1+edge):(M-edge);
phi_hw   = phi_hw(keep);  Sgphi_hw = Sgphi_hw(keep);  Sgi_hw = Sgi_hw(keep);  n_counts = n_counts(keep);
t        = t(keep);       M        = numel(t);
fprintf('Fs_dec = %.6f MHz, samples (trimmed) M = %d  (DECIM=%d)\n', Fs_dec/1e6, M, DECIM);

%% ---------------- build effective (decimated) kernel ----------------
% Full-rate kernel for one pulse
N_eff = floor(N_full/DECIM);                   % decimated samples per pulse
N_used = N_eff*DECIM;                          % truncate to a multiple of DECIM
tau  = (0:N_used-1)'/N_full;
g    = 1 ./ sqrt(tau + EPSG);  g(1) = 0;
dt_f = 1/N_full;
G    = sum(g)*dt_f;                            % ~2
K    = mu / G;
gdt  = g * dt_f;
h    = K * gdt;

% Block-sum by DECIM → effective decimated kernels
gdtD = sum(reshape(gdt, DECIM, N_eff), 1).';   % [N_eff x 1]
hD   = sum(reshape(h,   DECIM, N_eff), 1).';   % [N_eff x 1]
invSg = 1 / sum(gdtD);                         % decimated inv Sg

%% ---------------- pulse boundaries at decimated grid ----------------
pulses = detect_resets_from_Sgphi(Sgphi_hw, N_eff);
assert(~isempty(pulses), 'Could not detect pulse boundaries in Sg\phi.');
k0 = ceil(size(pulses,1)/2);
s0 = pulses(k0,1); e0 = pulses(k0,2);
midWin = (s0+round(0.1*N_eff)):(e0-round(0.1*N_eff));

clipEnd = pulses(end,2);        % last complete pulse end (decimated index)
clipIdx = 1:clipEnd;
phi_hw   = phi_hw(clipIdx);
Sgphi_hw = Sgphi_hw(clipIdx);
Sgi_hw   = Sgi_hw(clipIdx);
n_counts = n_counts(clipIdx);
t        = t(clipIdx);
M        = numel(t);

%% ---------------- per-sample ψ_true (decimated grid) ----------------
psi_true_n = zeros(M,1);
for r = 1:size(pulses,1)
    s = pulses(r,1); e = pulses(r,2);
    psi_true_n(s:e) = PSI_SEQ( 1 + mod(r-1, numel(PSI_SEQ)) );
end

%% ---------------- fit noise scale α on steady window ----------------
if FIT_NOISE_SCALE
    Ltmp = round(L_FIXED/DECIM);               % use decimated delay during fit
    fcost = @(a) rms_mid(phi_hw, sim_with_resets(N_eff,hD,gdtD, double(a)*double(n_counts), A, Ltmp, pulses, psi_true_n), midWin);
    alpha = fminsearch(fcost, 1.0, optimset('Display','off'));
else
    alpha = 1.0;
end

%% ---------------- choose delay L on decimated grid -----------------
switch lower(L_MODE)
    case 'fixed'
        L = round(L_FIXED/DECIM);
        phi_sm = sim_with_resets(N_eff,hD,gdtD, alpha*double(n_counts), A, L, pulses, psi_true_n);
    case 'grid'
        Ldec = unique(round(L_GRID/DECIM));
        best = inf; L = Ldec(1);
        for Lc = Ldec
            phi_c = sim_with_resets(N_eff,hD,gdtD, alpha*double(n_counts), A, Lc, pulses, psi_true_n);
            err   = rms_mid(phi_hw, phi_c, midWin);
            if err < best, best = err; L = Lc; phi_sm = phi_c; end
        end
    otherwise
        error('Unknown L_MODE');
end
fprintf('Delay (decimated grid): L=%d   (alpha=%.8g)\n', L, alpha);

%% ---------------- per-pulse metrics (same indices) ------------------
Sgphi_sm = sim_Sgphi_from_phi(gdtD,phi_sm,pulses);
Sgi_sm   = sim_Sgi_from_noise(gdtD,alpha*double(n_counts),A,phi_sm,pulses,psi_true_n);

[psi_e_hw, psi_i_hw, res_hw, tmid] = per_pulse(phi_hw, Sgphi_hw, Sgi_hw, pulses, invSg, A, Fs_dec);
[psi_e_sm, psi_i_sm, res_sm, ~   ] = per_pulse(phi_sm, Sgphi_sm, Sgi_sm, pulses, invSg, A, Fs_dec);

fprintf('Pulse @ t≈%.3f µs:\n', tmid(k0));
fprintf('  HW:  psi_end=%+.5f  psi_int=%+.5f  residual=%+.3e\n', psi_e_hw(k0), psi_i_hw(k0), res_hw(k0));
fprintf('  SIM: psi_end=%+.5f  psi_int=%+.5f  residual=%+.3e\n', psi_e_sm(k0), psi_i_sm(k0), res_sm(k0));

rms_abs = rms_mid(phi_hw, phi_sm, midWin);
rms_rel = rms_abs / max(1e-12, rms(phi_sm(midWin)));
fprintf('Steady window RMS error: %.3e rad  (relative %.3g)   [L=%d decimated]\n', rms_abs, rms_rel, L);

%% ---------------- plots ----------------
figure('Color','w','Position',[60 60 1200 640]);
plot(t*1e6, phi_hw,'b','LineWidth',1.05); hold on;
plot(t*1e6, phi_sm,'r--','LineWidth',1.0);
arrayfun(@(i) xline(t(pulses(i,1))*1e6,'k:'), 1:size(pulses,1));
plot(t*1e6, psi_true_n,'c:','LineWidth',1.0);
grid on; xlabel('Time [\mus]'); ylabel('\phi [rad]');
title(sprintf('\\phi: FPGA (blue) vs MATLAB (red)   [DECIM=%d, L=%d]', DECIM, L));
legend('FPGA','MATLAB sim','pulse start','\psi_{true}','Location','best');
end

%% ===================== helpers =====================
function [t, x, Fs] = load_trace(fname)
    T = readmatrix(fname); A = T(:,1:2);
    A = A(all(isfinite(A),2),:);
    t0 = A(:,1); x = A(:,2);
    dt = diff(t0); pos = dt(dt>0);
    meddt = median(pos); assert(~isempty(meddt) && isfinite(meddt),'Cannot infer dt');
    t  = t0; off=0; prev=t0(1);
    for i=2:numel(t0)
        if t0(i)<=prev, off = off + prev + meddt; end
        t(i) = t0(i)+off; prev=t0(i);
    end
    Fs = 1/meddt;             % this is the post-decimation Fs
end

function pulses = detect_resets_from_Sgphi(S, N_eff)
    dS = [0; diff(S)];
    thr0 = 0.02*max(abs(S)+eps);
    thrD = -5*median(abs(dS)+eps);
    idx  = find( (abs(S) < thr0) & (dS < thrD) );
    if isempty(idx)
        M = numel(S); b = 1:N_eff:M; e = min(b+N_eff-1,M);
        keep = (e-b+1)==N_eff; pulses=[b(keep)',e(keep)']; return;
    end
    min_gap = round(0.3*N_eff);
    idx = idx([true; diff(idx)>min_gap]);
    M = numel(S); S0 = idx(:); E0 = min(S0+N_eff-1, M);
    keep = (E0-S0+1)==N_eff; pulses = [S0(keep), E0(keep)];
end

function phi = sim_with_resets(N_eff,hD,gdtD,noise,A,L,pulses,psi_n)
    M = numel(noise); phi = zeros(M,1);
    k_idx = zeros(M,1);
    for r=1:size(pulses,1)
        s=pulses(r,1); e=pulses(r,2); k_idx(s:e) = (1:(e-s+1))';
    end
    dq = zeros(max(L,1),1); phi_cmd = 0;
    for n=1:M
        k  = k_idx(n); if k==0, k = 1 + mod(n-1,N_eff); end
        pa = dq(1);
        i  = A*sin(psi_n(n) - pa) + noise(n);
        phi_cmd = phi_cmd + hD(k) * i;     % decimated effective step
        if L>0, dq(1:end-1)=dq(2:end); dq(end)=phi_cmd; else, dq(1)=phi_cmd; end
        phi(n) = dq(1);
    end
end

function S = sim_Sgphi_from_phi(gdtD,phi,pulses)
    M = numel(phi); S = zeros(M,1);
    for r=1:size(pulses,1)
        s=pulses(r,1); e=pulses(r,2); k = (1:(e-s+1))';
        S(s:e) = cumsum(gdtD(k).*phi(s:e));
    end
end

function S = sim_Sgi_from_noise(gdtD,noise,A,phi,pulses,psi_n)
    M = numel(phi); S = zeros(M,1);
    for r=1:size(pulses,1)
        s=pulses(r,1); e=pulses(r,2); k = (1:(e-s+1))';
        i = A*sin(psi_n(s:e) - phi(s:e)) + noise(s:e);
        S(s:e) = cumsum(gdtD(k).*i);
    end
end

function [psi_end, psi_int, residual, tmid] = per_pulse(phi, sgphi, sgi, pulses, invSg, A, Fs)
    wrap = @(x) mod(x+pi,2*pi)-pi;
    P = size(pulses,1); psi_end=nan(P,1); psi_int=psi_end; residual=psi_end; tmid=psi_end;
    for r=1:P
        s=pulses(r,1); e=pulses(r,2);
        psi_end(r)  = wrap(phi(e));
        psi_int(r)  = wrap( (sgphi(e) + sgi(e)/A) * invSg );
        residual(r) = (sgi(e) * invSg) / A;
        tmid(r)     = ((s+e)/2 - 1)/Fs*1e6;
    end
end

function v = rms_mid(x, y, idx)
    v = rms(x(idx) - y(idx));
end
