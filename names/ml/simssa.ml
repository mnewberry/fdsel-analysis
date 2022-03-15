open Fdsel_lib ;;

let stof = float_of_string 
  and aofl = Array.of_list and atol = Array.to_list ;;
let simlabel = Sys.argv.(1)
let g_e = int_of_string (Sys.argv.(2)) ;;
let seed = int_of_string (Sys.argv.(3)) ;;
let centype = "CENSORED" ;;

let () = Rand.init seed ;;

let gtcs = parse_timeseries_gtc "inp/timeseries-ssaCd35dCE.tsv" ;;
let ts = gen_type_count_to_timeseries gtcs ;;
let true_uds = annupdate_data ts ;;

(* Simulate a wright-fisher expectation *) 
let sim_pops model indty =
  let wf_update (gen, true_ud) (maxty, ts) =
    (* Lift the true pop size and mutant counts from the real update data *)
    let mucs = Mu.filtmap (fun (_, (ic, fc)) ->
      if ic == 0 then Some fc else None) true_ud in
    let fps = pop_size_prime_a true_ud in
    let (maxty, mupop) = idpop_of_counts_string maxty mucs in
    (* Sample non-mutant population from last population *)
    let (last_gen, last_pop) = List.hd ts in
    let (lts, lcs) = Mu.unzip last_pop in
    let ips = Mu.sum lcs in 
    let iswt = map (fun (ty, ct) -> indty ips (ty, (ct, 0)) = 0) last_pop in
    let piis = List.tl (* ignore pimut *) (model last_pop) in
    let xis = atol (Gsl.Randist.multinomial Rand.rng fps (aofl piis)) in
    let pop = Mu.zip lts xis in
    let (maxty, pop) = Mu.rec_n (g_e - 1) (neutral_wf fps) (maxty, pop) in
    (* Introduce the true mutants and censor types. Note that this method of
       computing the censored counts preserves the number of wt individuals. *)
    let pop = Mu.filtmap2 (fun (ty, ct) wt ->
      if (not wt) && ct>=5 then Some (ty, ct) else None) pop iswt in
    let pop = pop @ mupop in
    let cenc = fps - pop_size_ts pop in
    let pop = (centype, cenc) :: pop in
    (maxty, (last_gen + 1, pop) :: ts) in
  (* get initial uncensored population from the data *)
  let (init_gen, init_uds) = List.hd true_uds in
  let init_pop = Mu.filtmap (fun (ty, (ct, _)) ->  match ct with
      0 -> None | _ -> Some (ty, ct)) init_uds in
  let (maxty, pops) = Mu.fold wf_update ("0", [(0, init_pop)]) true_uds in
  rev pops ;;

let nm_wf_exp wf pop =
  let (_, cts) = Mu.unzip pop in
  let ps = Mu.sum cts in
  let wis = map (fun ct -> wf (norm ct ps) *. norm ct ps) cts in
  let sum_wi = sumf wis in
  map (fun w -> (w /. sum_wi)) wis ;;

let ssaK () = 
  let indf = fun _ -> 0 in
  let swt = 0.0 in
  let pss = pop_sizes_auds_prime true_uds in
  let minfq = max_censored_freq 4 pss in
  let (nbins, indty) = parse_indty_args 1 indf (Some centype) (Some minfq) in
  let sels = aofl [swt;0.] in
  let model = wright_fisher_exp_ty 0. (w_pwcty indty sels) in
  sim_pops model indty ;;

let ssaFD () =
  let params = parse_params "inp/inf-ssaCd35dCE2Bfu.params.tsv" in
  let breaks = parse_nlsv stof "inp/inf-ssaCd35dCE2Bfu.breaks.tsv" in
  let minfreq = 0.00001 and swt = 0.0281510348703 in
  let indf = explicit_indf breaks in
  let (nbins, indty) = parse_indty_args (Array.length (fst params.sels))
    indf (Some centype) (Some minfreq) in
  let sels = aofl (swt :: (atol (fst params.sels))) in
  let model = wright_fisher_exp_ty 0. (w_pwcty indty sels) in
  sim_pops model indty ;;

let () = 
  let res = match simlabel with
      "FD" -> ssaFD ()
    | "K" -> ssaK () 
    | _ -> failwith "model must be FD or K" in
  fprint_pops_tsv Mu.identity
    (open_out 
      (Printf.sprintf "out/timeseries-ssa%s%dd35dCE.%03d.tsv" 
        simlabel g_e seed))
    "" res ;;
