//! Print every stage of the WOLF and PARSI paths for cd1 -> cd2, normalized for
//! comparison with DaliLite's serialcompare stage outputs (validation tooling).
//! Usage: stages <dat_dir> <cd1> <cd2>
use dali_core::pipeline::{dccp_to_dalicon_record, filter95_to_dalicon_records};
use dali_core::{dalicon, dp, filter95, parsi, wolf, AlignmentBlock, DccpEntry, ProteinStore};

fn blocks(b: &[AlignmentBlock]) -> String {
    b.iter().map(|x| format!("{}-{}:{}-{}", x.l1, x.r1, x.l2, x.r2)).collect::<Vec<_>>().join(" ")
}
fn dp_line(stage: &str, e: &DccpEntry) -> String {
    format!("{} {} {} {:.1} {:.1} : {}", stage, e.cd1, e.cd2, e.score, e.zscore, blocks(&e.blocks))
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (dir, cd1, cd2) = (&a[1], &a[2], &a[3]);
    let store = ProteinStore::new(dir);

    // WOLF path
    let d1 = wolf::setup_protein(format!("{}/{}.dat", dir, cd1)).unwrap();
    if d1.protein.nseg > 2 {
        let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
        wolf::load_protein(&mut grid, &d1);
        let d2 = wolf::setup_protein(format!("{}/{}.dat", dir, cd2)).unwrap();
        if let Some(w) = wolf::wolf_compare(&d1, &d2, &grid) {
            println!("wolf {} {} : {}", w.cd1, w.cd2, blocks(&w.blocks));
            let dpin = vec![(w.cd1.clone(), w.cd2.clone(), w.blocks.clone())];
            let dpo = dp::run_dp(&dpin, &store, dp::ZCUT_DEFAULT);
            for e in &dpo { println!("{}", dp_line("wolfdp", e)); }
            let recs: Vec<_> = dpo.iter().map(dccp_to_dalicon_record).collect();
            let dal = dalicon::run_dalicon(&recs, &store);
            for (x, y, _n, b) in &dal { println!("wolfdalicon {} {} : {}", x, y, blocks(b)); }
            let fin: Vec<_> = dal.into_iter().map(|(x, y, _n, b)| (x, y, b)).collect();
            for e in dp::run_dp(&fin, &store, dp::ZCUT_DEFAULT) { println!("{}", dp_line("wolffinal", &e)); }
        }
    }

    // PARSI path
    let mut cache = None;
    let hits = parsi::dowork_parsi(cd1, cd2, dir, None, true, &mut cache);
    for h in &hits {
        let r: Vec<String> = h.ranges_cd1.iter().chain(h.ranges_cd2.iter()).map(|p| format!("{} {}", p.0, p.1)).collect();
        println!("parsi {} {} {} {} : {}", h.cd1, h.cd2, h.idom, h.score, r.join(" "));
    }
    if hits.is_empty() { return; }
    let f95 = filter95::run_filter95(&hits, &store, Some(1.0));
    for e in &f95 {
        println!("filter95 {} {} {:.3} {} {} {} : {}", e.cd1, e.cd2, e.zscore, e.idom, e.score, e.nseg,
                 e.ranges.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" "));
    }
    let recs = filter95_to_dalicon_records(&f95, 1.0, 2);
    for r in &recs {
        println!("pipe96 {} {} {} : {}", r.cd1, r.cd2, r.nblock, r.values.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" "));
    }
    let dal = dalicon::run_dalicon(&recs, &store);
    for (x, y, _n, b) in &dal { println!("parsidalicon {} {} : {}", x, y, blocks(b)); }
    let fin: Vec<_> = dal.into_iter().map(|(x, y, _n, b)| (x, y, b)).collect();
    for e in dp::run_dp(&fin, &store, dp::ZCUT_DEFAULT) { println!("{}", dp_line("parsifinal", &e)); }
}
