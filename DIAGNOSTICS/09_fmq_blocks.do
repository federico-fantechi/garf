*! 09_fmq_blocks.do  --  forest match quality by block (paper: 5.9, Table 10)
*  ---------------------------------------------------------------------------
*  FMQ_i = (S_i - cbar_i) / (Smax_i - cbar_i), where
*    S_i    (match_cooccurrence)  : mean co-occurrence with i's own donors
*    cbar_i (random_cooccurrence) : mean co-occurrence with random untreated units
*    Smax_i (best_cooccurrence)   : co-occurrence with the closest untreated unit
*  0 = donors no closer than random controls; 1 = every donor as close as the
*  best one the data contain. fmq_twin = (S_i - cbar_i)/(1 - cbar_i) is the
*  beta 1.2-1.4 definition (comparable only within the same N and forest).
*  A LOW FMQ has two readings:
*    high cbar_i            -> a typical unit, any control is about as close
*    low cbar_i, low S_i,
*    low overlap mass       -> a genuinely isolated unit
*  Always read FMQ together with ESS (07_ess_reuse.do).
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
import delimited "garf_output/blocks.csv", clear varnames(1)

tabstat fmq fmq_twin ess match_cooccurrence random_cooccurrence best_cooccurrence, ///
    stats(mean p10 p50 p90 min max) format(%5.3f)
count if missing(fmq) & matched == 1     // blocks with undefined FMQ (Smax = cbar)

* the five treated units with the lowest FMQ, with everything needed to read them
sort fmq
list treated_id fmq match_cooccurrence random_cooccurrence best_cooccurrence ess ///
     n_donors best_proximity_minmax overlap_mass in 1/5, noobs
