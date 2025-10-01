use qfilter::merkle;
use rand::{rng, Rng};
use std::collections::HashSet;

type Run = (u64, u64, u64, Vec<u8>);
fn run_contains(f: qfilter::Filter, hash: u64, run: Run) -> bool {
    let (q_bucket_idx, run_start_idx, run_end_idx, bytes) = run;
    let blocks_count = bytes.len() / f.block_byte_size();

    let mut is_contains = false;
    for block_num in 0..blocks_count {
        let block_start_bit = block_num as usize * f.block_byte_size();

        let block_bytes;
        if block_num > 0 {
            block_bytes = &bytes[block_start_bit..][..f.block_byte_size()];
        } else {
            block_bytes = &bytes[..f.block_byte_size()];
        }

        let reminder_size = f.rbits.get() as usize;
        let reminders = &block_bytes[1 + 8 + 8..][..8 * reminder_size];

        let mut run_start_intra_block_idx = 0;
        let mut run_end_intra_block_idx = 63;
        // single block
        if block_num == 0 && blocks_count == 1 {
            run_start_intra_block_idx = run_start_idx % 64;
            run_end_intra_block_idx = run_end_idx % 64;
            // first block in block sequence
        } else if block_num == 0 {
            run_start_intra_block_idx = run_start_idx % 64;
        // last block in block sequence
        } else if block_num == blocks_count - 1 {
            run_end_intra_block_idx = run_end_idx % 64;
        }
        for bucket_idx in run_start_intra_block_idx..=run_end_intra_block_idx {
            let bucket_start_bit = bucket_idx as u64 * f.rbits.get() as u64;
            let bucket_end_bit = bucket_start_bit + f.rbits.get() as u64;

            let start_u64 = (bucket_start_bit / 64) as usize;
            let num_rem_parts = 1 + (bucket_end_bit > ((start_u64 + 1) * 64) as u64) as usize;

            let extra_low = bucket_start_bit as usize - start_u64 * 64;
            let extra_high =
                (((start_u64 + 1) * 64).saturating_sub(bucket_end_bit as usize)) as usize;

            let rem_parts_bytes = &reminders[start_u64 * 8..][..num_rem_parts * 8];
            let rem_part = u64::from_le_bytes(rem_parts_bytes[..8].try_into().unwrap());
            let mut remainder = (rem_part << extra_high) >> (extra_high + extra_low);

            if let Some(rem_part) = rem_parts_bytes.get(8..16) {
                let remaining_bits = bucket_end_bit - ((start_u64 + 1) * 64) as u64;
                let rem_part = u64::from_le_bytes(rem_part.try_into().unwrap());
                remainder |= (rem_part & !(u64::MAX << remaining_bits))
                    << (f.rbits.get() as u64 - remaining_bits);
            }

            let quotient = q_bucket_idx << f.rbits.get();

            let run_hash = quotient | remainder;

            if run_hash == hash {
                is_contains = true;
            }
        }
    }

    is_contains
}

fn main() {
    let items_len = 100000;
    let mut f = qfilter::Filter::new(items_len, 0.005).unwrap();

    // Generate test data
    let rem_mask = (1 << f.rbits.get() as usize) - 1;
    let quot_mask = ((1 << f.qbits.get() as usize) - 1) << f.rbits.get() as usize;

    let collision_quot = rng().random::<u64>() & quot_mask;
    let mut items: Vec<u64> = vec![];
    for i in 0..items_len {
        if (i % 2 == 0 || i % 3 == 0) && i % 6 != 0 {
            let rem = rng().random::<u64>();
            let rem = rem & rem_mask;
            let item = collision_quot | rem;
            items.push(item);
        } else {
            let quot = rng().random::<u64>() & quot_mask;
            let rem = rng().random::<u64>() & rem_mask;
            let item = quot | rem;
            items.push(item);
        }
    }

    for &i in items.iter() {
        f.insert(i);
    }

    let unique_items = items.clone().into_iter().collect::<HashSet<u64>>();

    println!("unique items len: {:?}", unique_items.len());

    println!("quotient filter: {:?} (len: {}))", f, f.len());

    let mut found_items = 0;

    for (q_bucket_idx, run_start_idx, run_end_idx, bytes) in f.run_blocks() {
        let blocks_count = bytes.len() / f.block_byte_size();

        for block_num in 0..blocks_count {
            let block_start_bit = block_num as usize * f.block_byte_size();

            let block_bytes;
            if block_num > 0 {
                block_bytes = &bytes[block_start_bit..][..f.block_byte_size()];
            } else {
                block_bytes = &bytes[..f.block_byte_size()];
            }

            let reminder_size = f.rbits.get() as usize;
            let reminders = &block_bytes[1 + 8 + 8..][..8 * reminder_size];

            let mut run_start_intra_block_idx = 0;
            let mut run_end_intra_block_idx = 63;
            // single block
            if block_num == 0 && blocks_count == 1 {
                run_start_intra_block_idx = run_start_idx % 64;
                run_end_intra_block_idx = run_end_idx % 64;
                // first block in block sequence
            } else if block_num == 0 {
                run_start_intra_block_idx = run_start_idx % 64;
            // last block in block sequence
            } else if block_num == blocks_count - 1 {
                run_end_intra_block_idx = run_end_idx % 64;
            }
            for bucket_idx in run_start_intra_block_idx..=run_end_intra_block_idx {
                let bucket_start_bit = bucket_idx as u64 * f.rbits.get() as u64;
                let bucket_end_bit = bucket_start_bit + f.rbits.get() as u64;

                let start_u64 = (bucket_start_bit / 64) as usize;
                let num_rem_parts = 1 + (bucket_end_bit > ((start_u64 + 1) * 64) as u64) as usize;

                let extra_low = bucket_start_bit as usize - start_u64 * 64;
                let extra_high =
                    (((start_u64 + 1) * 64).saturating_sub(bucket_end_bit as usize)) as usize;

                let rem_parts_bytes = &reminders[start_u64 * 8..][..num_rem_parts * 8];
                let rem_part = u64::from_le_bytes(rem_parts_bytes[..8].try_into().unwrap());
                let mut remainder = (rem_part << extra_high) >> (extra_high + extra_low);

                if let Some(rem_part) = rem_parts_bytes.get(8..16) {
                    let remaining_bits = bucket_end_bit - ((start_u64 + 1) * 64) as u64;
                    let rem_part = u64::from_le_bytes(rem_part.try_into().unwrap());
                    remainder |= (rem_part & !(u64::MAX << remaining_bits))
                        << (f.rbits.get() as u64 - remaining_bits);
                }

                let quotient = q_bucket_idx << f.rbits.get();

                let hash = quotient | remainder;
                if items.contains(&hash) {
                    found_items += 1;
                } else {
                }
            }
        }
    }
    assert_eq!(unique_items.len(), found_items);

    let merkle_leaves = merkle::leaf_hashes_from_iter(f.run_blocks().map(|(_, _, _, bytes)| bytes));
    let merkle_root = merkle::compute_root(merkle_leaves.clone());

    let item_index = rng().random_range(0..items_len);
    if let Some(item) = items.get(item_index as usize) {
        let (item_q, item_r) = f.calc_qr(*item);
        println!("item quotient: {} ({:064b})", item_q, item_q);
        if let Some((index, item_potential_run)) =
            f.run_blocks().enumerate().find(|(i, run)| run.0 == item_q)
        {
            println!("item potential run {:?}", item_potential_run.0);
            if run_contains(f, *item, item_potential_run) == true {
                let item_membership_proof = merkle::build_proof_path(&merkle_leaves, index);
            };
        };
    };
}
