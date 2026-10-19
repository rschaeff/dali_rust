//! Print every stage of the SOAP path for cd1 -> each cd2, normalized for
//! comparison with DaliLite's serialcompare stage outputs (validation tooling).
//! Usage: soapstages <dat_dir> <cd1> <cd2> [<cd2> ...]
use dali_core::pipeline::dccp_to_dalicon_record;
use dali_core::{dalicon, dp, soap, AlignmentBlock, DccpEntry, ProteinStore};

fn blocks(b: &[AlignmentBlock]) -> String {
    b.iter().map(|x| format!("{}-{}:{}-{}", x.l1, x.r1, x.l2, x.r2)).collect::<Vec<_>>().join(" ")
}
fn dp_line(stage: &str, e: &DccpEntry) -> String {
    format!("{} {} {} {:.1} {:.1} : {}", stage, e.cd1, e.cd2, e.score, e.zscore, blocks(&e.blocks))
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (dir, cd1) = (&a[1], &a[2]);
    let store = ProteinStore::new(dir);
    let mut q = soap::SoapQuery::new(&store.get_protein(cd1).unwrap());
    let mut w = Vec::new();
    for cd2 in &a[3..] {
        let t = store.get_protein(cd2).unwrap();
        if let Some(r) = soap::soap_compare(&mut q, &t) {
            println!("soap {} {} : {}", r.cd1, r.cd2, blocks(&r.blocks));
            w.push((r.cd1, r.cd2, r.blocks));
        }
    }
    let dpo = dp::run_dp(&w, &store, dp::ZCUT_DEFAULT);
    for e in &dpo { println!("{}", dp_line("soapdp", e)); }
    let recs: Vec<_> = dpo.iter().map(dccp_to_dalicon_record).collect();
    let dal = dalicon::run_dalicon(&recs, &store);
    for (x, y, _n, b) in &dal { println!("soapdalicon {} {} : {}", x, y, blocks(b)); }
    let fin: Vec<_> = dal.into_iter().map(|(x, y, _n, b)| (x, y, b)).collect();
    for e in dp::run_dp(&fin, &store, dp::ZCUT_DEFAULT) { println!("{}", dp_line("soapfinal", &e)); }
}
