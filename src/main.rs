use qfilter::merkle;
use qfilter::{packed, run};

use std::io;

use itertools::Itertools;
use rand::{rand_core::block, rng, Rng};
use std::hash::Hash;

const CELESTIA_SAMPLE_SIZE: usize = 482;

#[derive(Debug)]
struct RunSpan {
    idx: u64,
    start: u64,
    end: u64,
}

fn main() {
    let items_len = 1000;
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

    println!("Packing blocks in 512 chunks");

    let mut packed_qf: Vec<packed::PackedBlocks> = vec![];

    // let sample_block_num = CELESTIA_SAMPLE_SIZE / f.block_byte_size();
    let sample_block_num = 1;

    let total_blocks = f.total_blocks().get() as u64;
    for block_range_start in (0..2).step_by(sample_block_num) {
        let block_range_end = (block_range_start + sample_block_num as u64).min(total_blocks);

        let mut packed_runs: Vec<run::Run> = vec![];
        let mut packed_blocks = packed::PackedBlocks::default();
        packed_blocks.block_offset = block_range_start;
        packed_blocks.blocks_num = (block_range_end - block_range_start) as u8;
        packed_blocks.rbits = f.rbits.get() as u8;
        packed_blocks.qbits = f.qbits.get() as u8;

        for block_num in block_range_start..block_range_end {
            println!(
                "________________________ block {} ------------------------",
                block_num
            );
            println!("packed blocks: {:?}", packed_blocks);
            let mut buffer = String::new();
            let (mut runs, shifted) = f.block_contains_runs(block_num);
            let runs_idx: Vec<u64> = runs.iter().map(|run| run.q_bucket_idx).collect();
            println!("runs idx: {:?}", runs_idx);

            packed_runs.append(&mut runs);
            let block = f.block(block_num);
            let block_bytes: Vec<u8> = f.block_bytes_with_r(block_num).to_vec();

            println!("offset: {}", block.offset);
            println!("occupiends [{:b}]", block.occupieds);
            println!("runends [{:b}]", block.runends);

            println!(
                "block {block_num} contains {:?} [{:?}, {:?}]",
                runs.len(),
                block_num * 64,
                block_num * 64 + 64,
            );

            if block_num == block_range_end {
                packed_blocks.shifted = shifted;
            }
            packed_blocks.first_runstart_idx = runs_idx
                .first()
                .map(|&current_block_first_runstart_idx| {
                    packed_blocks
                        .first_runstart_idx
                        .map(|first_runstart_idx| {
                            current_block_first_runstart_idx.min(first_runstart_idx)
                        })
                        .unwrap_or(current_block_first_runstart_idx)
                })
                .or(packed_blocks.first_runstart_idx);

            packed_blocks.buffer.append(&mut block_bytes.to_vec());
            let from_packed_block = packed_blocks.raw_block(block_num);
            println!("offset: {}", from_packed_block.offset);
            println!("occupiends [{:b}]", from_packed_block.occupieds);
            println!("runends [{:b}]", from_packed_block.runends);

            println!(
                "________________________ block {} end --------------------",
                block_num
            );
        }
        let mut contains_count = 0;
        for item in items.iter() {
            println!("^^^^^^^^^^^^^^^^^^^^^^^^^^^");
            if packed_blocks.contains(item) {
                contains_count += 1;
                println!("packed blocks contains {} items", contains_count);
            }
        }

        let run_spans: Vec<RunSpan> = packed_runs
            .iter()
            .filter_map(|run| {
                if let Some((start, end)) = run.start_idx.zip(run.end_idx) {
                    Some(RunSpan {
                        idx: run.q_bucket_idx,
                        start,
                        end,
                    })
                } else {
                    None
                }
            })
            .collect();

        let count_buckets: u64 = run_spans
            .iter()
            .map(|run_span| {
                println!("count bucket: [{}, {}]", run_span.start, run_span.end);
                let buckets_num = run_span.end - run_span.start + 1;
                buckets_num
            })
            .sum();
        let count_runs = packed_runs.len();
        println!("count packed runs : {}", count_runs);

        println!("packed blocks contains {} items", contains_count);
        println!("count packed runs bukcets: {}", count_buckets);

        // packed_runs.iter().for_each(|run| {
        //     // let runs_span = packed_blocks.runs_span();
        //
        //     // if run.q_bucket_idx > runs_span.0 && run.q_bucket_idx < runs_span.1 {
        //     //     println!("run end for {}", run.q_bucket_idx);
        //     //     println!("run end: {}", packed_blocks.run_end(run.q_bucket_idx));
        //     // }
        //
        // });
        packed_qf.push(packed_blocks);
    }

    println!("Packed blocks in 512 chunks");
}
