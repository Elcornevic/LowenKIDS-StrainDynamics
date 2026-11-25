# LowenKIDS-StrainDynamics
Bioinformatic and phylogenetic pipeline (R, Python, ETE 3) analyzing strain-level genetic distances, longitudinal persistence/turnover, and household transmission across 932 gut microbiome SGBs in the LöwenKIDS cohort using StrainPhlAn4 outputs.

```mermaid
flowchart TD
    %% ==================== TIER 1: INGESTION, QC & METADATA ====================
    subgraph T1["Tier 1 — Data Ingestion, QC & Master Matrix (932 SGBs)"]
        direction LR
        M1["M1: Raw Inputs (Bash)<br/>932 SGB dirs (t__SGB*_tree/)<br/>RAxML_bestTree.tre, .aln, .mutation<br/>.polymorphic, .info, _mutation_rates<br/>LöwenKIDS metadata & CHOCOPhlAn tax"]
        M2["M2: Count & .aln Audit (R/Bash)<br/>Exact *concatenated.aln (no .reduced)<br/>Strip '>' & '_filtered' sample IDs<br/>Reconcile .polymorphic, .mutation, .aln<br/>Parse marker count & L_aln"]
        M3["M3: Taxonomy & Metadata (R)<br/>Strip 't__' & '_group' fallback<br/>Join sample_id, household_id,<br/>timestamp & role ('Kid')<br/>Compute n_households & n_kids_multi_hh"]
        M4["M4: Power Tiers (R)<br/>Tier A: All 932 SGBs<br/>Tier B: n_aln >= 10, HH >= 3<br/>Tier C: n_kids_multi_hh >= 3, Δt > 0<br/>Flag replicates (Δt = 0)"]
        M1 --> M2 --> M3 --> M4
    end

    %% ==================== TIER 2: DISTANCES & ETE 3 ====================
    subgraph T2["Tier 2 — Pairwise Distances, Polymorphism & ETE 3 Phylogenetics"]
        direction LR
        M5["M5: Mutation Pairs (Python/pandas)<br/>Melt *.mutation -> n(n-1)/2 pairs<br/>Compute d_mut(i, j)<br/>Cross-check .aln SNP distance d_aln"]
        M6["M6: Polymorphism QC (Python/R)<br/>Parse percentage_of_polymorphic_sites<br/>& marker stats from *.polymorphic<br/>Flag mixed-strain / high-P_poly samples"]
        M7["M7: ETE 3 Patristic (Python/ete3)<br/>Load RAxML_bestTree.tre<br/>Midpoint root & strip '_filtered'<br/>Compute d_pat(i, j) & ρ(d_mut, d_pat)"]
        M8["M8: ETE 3 Topology (Python/ete3)<br/>Terminal branch lengths l_term<br/>Nearest-neighbor HH match NN_HH<br/>Monophyly & clade clustering per HH"]
        M5 --> M6 --> M7 --> M8
    end

    %% ==================== TIER 3: PAIRS & CUTOFF CALIBRATION ====================
    subgraph T3["Tier 3 — Relational Pair Taxonomy & Cutoff Calibration"]
        direction LR
        M9["M9: Relational Axes (Python/R)<br/>Within-HH vs Between-HH<br/>Same time (Δt=0) vs Longitudinal (Δt>0)<br/>Role pairs (Kid-Kid, Kid-HH, Inter-HH)"]
        M10["M10: Cutoff τ Sweep (R/Python)<br/>Calibrate τ_mut & τ_pat from W_s vs B_s<br/>Indicator S_ij(τ) (d_mut <= τ & d_pat <= τ)<br/>Sensitivity sweep across τ grid"]
        M11["M11: Mutation Rates (R/Python)<br/>Expected drift E[d(Δt)] = μ_s * Δt + d_0<br/>Excess divergence Δd_excess<br/>Separate drift vs strain turnover"]
        M9 --> M10 --> M11
    end

    %% ==================== TIER 4: QUESTION 1 (LONGITUDINAL) ====================
    subgraph T4["Tier 4 — Q1: How Strains Stay in One Subject Over Time or Change?"]
        direction LR
        Q1_1["Q1.1: Transitions (Python/R)<br/>Ordered series m_1 -> ... -> m_k<br/>Persistence: d(m_r, m_{r+1}) <= τ<br/>Replacement: d(m_r, m_{r+1}) > τ<br/>PR_s & TR_s = 1 - PR_s"]
        Q1_2["Q1.2: Resurgence (Python/R)<br/>Absent at t_mid in .aln, present at t_3<br/>Resurgence: d(m_{t1}, m_{t3}) <= τ<br/>Re-colonization: d(m_{t1}, m_{t3}) > τ<br/>Resurgence fraction R_res,s"]
        Q1_3["Q1.3: Time & Survival (R)<br/>d_ij ~ β_0 + β_1 Δt_ij + (1|HH) + (1|SGB)<br/>Kaplan-Meier strain retention curves<br/>Strain half-life t_1/2 per SGB"]
        Q1_4["Q1.4: In-Host Evolution (Python/R)<br/>SNV accumulation r_inhost = d_mut / Δt<br/>Polymorphism P_poly before turnover<br/>Co-colonization vs replacement"]
        Q1_1 --> Q1_2 --> Q1_3 --> Q1_4
    end

    %% ==================== TIER 5: QUESTION 2 (HOUSEHOLDS) ====================
    subgraph T5["Tier 5 — Q2: How Similar or Different Are Strains Between Households?"]
        direction LR
        Q2_1["Q2.1: W_s vs B_s Tests (R/scipy)<br/>Mann-Whitney U & 9,999 permutations<br/>BH-FDR adjusted q-values<br/>Δ_med, R_med & Cliff's delta δ"]
        Q2_2["Q2.2: Shared Strains (R/Python)<br/>P_intra(s) vs P_inter(s) under τ<br/>HH sharing odds ratio OR_HH(s)<br/>Δt-matched household comparisons"]
        Q2_3["Q2.3: Role Network (R/igraph)<br/>Kid-Parent, Kid-Sibling, Adult-Adult<br/>Vertical vs horizontal sharing<br/>HH size & Kid age effects"]
        Q2_4["Q2.4: ETE 3 HH Clades (ete3)<br/>W_s vs B_s on patristic d_pat<br/>Sub-threshold / monophyletic HH clades<br/>Phylogenetic HH signal"]
        Q2_1 --> Q2_2 --> Q2_3 --> Q2_4
    end

    %% ==================== TIER 6: SYNTHESIS & DELIVERABLES ====================
    subgraph T6["Tier 6 — Cross-SGB Synthesis (932 SGBs) & Thesis Deliverables"]
        direction LR
        S1["S1: Master 932 Matrix<br/>Coverage, Q1, Q2, ETE 3,<br/>mutation rates & P_poly"]
        S2["S2: 4-Quadrant Typology<br/>(PR_s vs OR_HH) ecological<br/>profiles across SGBs"]
        S3["S3: Taxa & Traits<br/>Phylum/Family/Genus tests<br/>& genomic predictors"]
        S4["S4: Robustness & Outputs<br/>τ & Tier sensitivity, Git repo,<br/>figures & HZI deck"]
        S1 --> S2 --> S3 --> S4
    end

    %% ==================== CONNECTIONS BETWEEN TIERS ====================
    M4 --> M5
    M4 --> M7
    M8 --> M9
    M6 --> M9
    M11 --> Q1_1
    M11 --> Q2_1
    Q1_4 --> S1
    Q2_4 --> S1
```
