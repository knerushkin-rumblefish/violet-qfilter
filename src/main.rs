use qfilter::merkle;

use rand::{rng, Rng};
use std::hash::Hash;

fn main() {
    let items_len = 100000;
    let mut f = qfilter::Filter::new(items_len, 0.005).unwrap();

    let mut items: Vec<[u8; 32]> = vec![];
    for _ in 0..items_len {
        let item = rng().random::<[u8; 32]>();
        items.push(item);
    }

    for &i in items.iter() {
        if let Err(e) = f.insert(i) {
            eprintln!("item insertion failed: {e}");
        }
    }

    let merkle_leaves = merkle::leaf_hashes_from_iter(f.run_blocks().map(|run| run.buffer));
    let merkle_root = merkle::compute_root(merkle_leaves.clone());

    let item_index = rng().random_range(0..items_len);
    if let Some(item) = items.get(item_index as usize) {
        let (item_q_bucket_idx, _) = f.calc_item_qr(*item);
        println!(
            "item quotient: {} ({:064b})",
            item_q_bucket_idx, item_q_bucket_idx
        );
        println!("run count: {}", f.run_blocks().count());
        if let Some((index, item_potential_run)) = f.run_blocks().enumerate().find(|(_, run)| {
            //println!("run {}", run.q_bucket_idx);
            run.q_bucket_idx == item_q_bucket_idx
        }) {
            println!("item potential run {:?}", item_potential_run.q_bucket_idx);
            if f.contains(item) {
                println!("item is contains in filter")
            }
            if item_potential_run.contains(item) {
                println!("item is contains in run");
                let item_membership_proof = merkle::build_proof_path(&merkle_leaves, index);
            };
        };
    };
}
