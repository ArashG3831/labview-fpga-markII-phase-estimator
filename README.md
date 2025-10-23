# Mark-II Weighted Phase Estimator (Real-Time, 125 MHz)

**What**: Per-pulse weighted estimator on PXIe-7972R + NI-5782R. Weights (two vectors, **N=6250** each) are precomputed in MATLAB, quantized to **Q4.14**, shipped via two **Host→Target DMAs** into dual-clock BRAM. A **Weights Ready** gate releases the 125 MHz SCTL. Achieves steady-window parity ≈ **2.73×10⁻³ rad RMS** vs MATLAB. :contentReference[oaicite:7]{index=7}

**Why precompute weights?** Division and √ on FPGA would bloat latency; we treat the card as a high-rate MAC engine with fixed coeffs. :contentReference[oaicite:8]{index=8}

**Hardware/Tools**: Same as the IIR repo (PXIe-7972R, NI-5782R, LabVIEW FPGA, MATLAB). :contentReference[oaicite:9]{index=9}

**Quickstart**
1) Place the two CSVs in `/data` (**6250 lines each**, no headers). :contentReference[oaicite:10]{index=10}  
2) Run host VI: read CSV → quantize to **Q4.14 (18-bit word, 4 int)** → write both **Host→Target DMA** FIFOs. Timeout ~5 s. :contentReference[oaicite:11]{index=11}  
3) FPGA loader (fabric clock) writes BRAM, asserts **Weights Ready**, SCTL starts; streams per-pulse summaries & timeseries to host. :contentReference[oaicite:12]{index=12}  
4) Run `matlab/validate_markII.m` for overlays (HW vs MATLAB) and RMS parity. Expect ≈ **2.7e-3 rad** on a steady window. :contentReference[oaicite:13]{index=13}

**Notes**
- Signals carried as **Q2.13**; weights as **Q4.14**; products widened on DSP P path; single cast near outputs. :contentReference[oaicite:14]{index=14}
- If you decimate streams, set host waveform **dt = 8 ns × Ndec**. :contentReference[oaicite:15]{index=15}

**Reference**: The full UNSW report details the DMA names, BRAM config, loader handshake, and AO drive path. :contentReference[oaicite:16]{index=16}
