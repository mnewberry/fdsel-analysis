open Fdsel_lib

let wright_fisher_exp mu wf pop =
  let pimut = mu in
  let ps = pop_size_ts pop in
  let wis = map (fun (ty, ct) -> wf ps (ty, ct) *. norm ct ps) pop in
  let sum_wi = sumf wis in
  let piis = map (fun w -> (w /. sum_wi) *. (1. -. pimut)) wis in
  pimut :: piis

let wright_fisher time_rescaling mu wf ps (maxty, pop) =
  let mult_piis = aofl (wright_fisher_exp mu wf pop) in
  let xps = atol (Gsl.Randist.multinomial Rand.rng ps mult_piis) in
  let mut_ct, nm_counts = List.hd xps, List.tl xps in
  let mut_counts = Array.make mut_ct 1 in
  let maxty, mut_pop = idpop_of_counts maxty (atol mut_counts) in
  let result = (maxty,
    Mu.fold2 (fun (ty,_) ct pop -> (ty,ct)::pop) mut_pop pop nm_counts) in
  Mu.rec_n (time_rescaling - 1) (neutral_wf ps) result

let wright_fisher_ time_rescaling lambda mu wf ps (maxty, pop) =
  let mult_piis = aofl (wright_fisher_exp mu wf pop) in
  let xps = atol (Gsl.Randist.multinomial Rand.rng ps mult_piis) in
  let mut_tot, nm_counts = List.hd xps, List.tl xps in
  let mut_ntypes = if mut_tot > 1
    then 1 + (Gsl.Randist.binomial Rand.rng lambda (mut_tot - 1))
    else mut_tot in
  let mut_counts = Array.map ((+) 1)
    (Gsl.Randist.multinomial Rand.rng (mut_tot - mut_ntypes)
    (Array.make mut_ntypes (1. /. fl mut_ntypes))) in
  let maxty, mut_pop = idpop_of_counts maxty (atol mut_counts) in
  let result = (maxty,
    Mu.fold2 (fun (ty,_) ct pop -> (ty,ct)::pop) mut_pop pop nm_counts) in
  Mu.rec_n (time_rescaling - 1) (neutral_wf ps) result

let ps = 20000 and mu = 0.0003 and ngens = 50000
let monopop = monomorphic_idpop ps

let simulate wf filename =
  let outch = open_out filename in
  let model = wright_fisher_ 1 0.5 mu wf in
  let init = simulate ~seed:(Some 1) ~demog:(Fixed(ps, ps))
    agg_none (snd monopop) model monopop in
  let init = let (mint, maxt) = Mu.min_max (map fst init) in
    (maxt - mint + 1, map (fun (ty,ct) -> (ty - mint + 1, ct)) init) in
  fprint_pops_tsv_header outch "" ;
  ignore (simulate ~seed:(Some 0) ~demog:(Fixed (ps, ngens))
    (agg_all_print (fprint_pops_tsv_inc string_of_int outch "") 1)
    (agg_all_knil (snd init)) model init)
   
let fpar ps (ty, ct) = let x = norm ct ps in
  1. +. 0.01 *. (-4. -. log10 x) *. (log10 x +. 2.)

let flin ps (ty, ct) = let x = norm ct ps in
  1.012 -. 0.006 *. (log10 x +. 4.)

let fparlin ps (ty, ct) =
  if ty mod 3 = 0 then fpar ps (ty, ct) else flin ps (ty, ct)

let flinoff ps (ty, ct) =
  if ty mod 3 = 0 then flin ps (ty, ct) +. 0.005 else flin ps (ty, ct)

let () =
  simulate fpar "timeseries-simpar.tsv" ;
  simulate flin "timeseries-simlin.tsv" ;
  simulate fparlin "timeseries-simparlin.tsv" ;
  simulate flinoff "timeseries-simlinoff.tsv"
