use crate::{BitExt, Block as BytesBlock, Filter, Membership};

use std::hash::Hash;

const BLOCK_OFFSET_SIZE: usize = 1;
const BLOCK_OCCUPIEDS_SIZE: usize = 8;
const BLOCK_RUNENDS_SIZE: usize = 8;

const BLOCK_REMINDERS_SLOTS_SIZE: usize = 64;

#[derive(Default, Debug, Clone)]
pub struct Block {
    pub first_q_bucket_idx: u64,
    pub first_q_bucket_runstart_idx: Option<u64>,

    pub buffer: Vec<u8>,

    pub rbits: u8,
    pub qbits: u8,

    pub next: Option<Box<Block>>,
}

impl Block {
    #[inline]
    fn block_header_size(&self) -> usize {
        BLOCK_OFFSET_SIZE + BLOCK_OCCUPIEDS_SIZE + BLOCK_RUNENDS_SIZE
    }

    #[inline]
    pub fn block_byte_size(&self) -> usize {
        self.block_header_size() + BLOCK_REMINDERS_SLOTS_SIZE * self.rbits as usize / 8
    }

    #[inline]
    fn block_bytes(&self) -> &[u8] {
        &self.buffer[..self.block_header_size()]
    }

    #[inline]
    fn total_buckets(&self) -> u64 {
        1 << self.qbits
    }

    // TODO: this is part of QF logic so should be in appropriate trait
    #[inline]
    fn calc_qr(&self, hash: u64) -> (u64, u64) {
        let q = (hash >> self.rbits) & ((1 << self.qbits) - 1);
        let r = hash & ((1 << self.rbits) - 1);
        (q, r)
    }

    #[inline]
    fn raw_block(&self) -> BytesBlock {
        let block_bytes: &[u8; 1 + 8 + 8] = &self.block_bytes().try_into().unwrap();

        BytesBlock {
            offset: block_bytes[0] as u64,
            occupieds: u64::from_le_bytes(
                block_bytes[BLOCK_OFFSET_SIZE..BLOCK_OFFSET_SIZE + BLOCK_OCCUPIEDS_SIZE]
                    .try_into()
                    .unwrap(),
            ),
            runends: u64::from_le_bytes(
                block_bytes[BLOCK_OFFSET_SIZE + BLOCK_OCCUPIEDS_SIZE
                    ..BLOCK_OFFSET_SIZE + BLOCK_OCCUPIEDS_SIZE + BLOCK_RUNENDS_SIZE]
                    .try_into()
                    .unwrap(),
            ),
        }
    }

    pub fn run_end(&self, hash_bucket_idx: u64) -> Option<u64> {
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();
        assert!(self.contains_q_bucket_idx(hash_bucket_idx));

        let hash_bucket_block_idx = hash_bucket_idx / 64;
        assert_eq!(hash_bucket_block_idx, self.first_q_bucket_idx / 64);

        let hash_bucket_intrablock_offset = hash_bucket_idx % 64;

        let current_block = self.raw_block();
        let hash_bucket_intrablock_rank = current_block
            .occupieds
            .popcnt(..=hash_bucket_intrablock_offset);
        // No occupied buckets all the way to bucket_intrablock_offset
        // which also means hash_bucket_idx isn't occupied
        if hash_bucket_intrablock_rank == 0 {
            // TODO: WTF???
            return if current_block.offset <= hash_bucket_intrablock_offset {
                // hash_bucket_idx points to an empty bucket unaffected by block offset,
                // thus end == start
                Some(hash_bucket_idx)
            } else {
                // hash_bucket_idx fall within the section occupied by the offset,
                // thus end == last bucket of offset section
                // TODO: WTF???
                Some((hash_bucket_block_idx * 64 + current_block.offset - 1) % self.total_buckets())
            };
        }

        // Must search runends to figure out the end of the run
        // Is this block or next one?
        // if next we should follow to extension
        let mut runend_block_idx = hash_bucket_block_idx + current_block.offset / 64;
        assert_eq!(runend_block_idx, self.first_q_bucket_idx / 64);

        let mut runend_ignore_bits = current_block.offset % 64;
        let mut runend_block = current_block;
        let mut runend_rank = hash_bucket_intrablock_rank - 1;

        let mut runend_block_offset = runend_block
            .runends
            .select(runend_ignore_bits.., runend_rank);

        if let Some(runend_block_offset) = runend_block_offset {
            let runend_idx = runend_block_idx * 64 + runend_block_offset;
            return Some(runend_idx.max(hash_bucket_idx) % self.total_buckets());
        }

        loop {
            if let Some(next_block) = self.next.clone() {
                runend_rank -= runend_block.runends.popcnt(runend_ignore_bits..);
                runend_block_idx += 1;
                runend_ignore_bits = 0;
                runend_block = next_block.raw_block();
                runend_block_offset = runend_block
                    .runends
                    .select(runend_ignore_bits.., runend_rank);

                if let Some(runend_block_offset) = runend_block_offset {
                    let runend_idx = runend_block_idx * 64 + runend_block_offset;
                    return Some(runend_idx.max(hash_bucket_idx) % self.total_buckets());
                }
            }
        }
    }

    //
    // #[inline]
    // fn run_start(&self, hash_bucket_idx: u64) -> u64 {
    //     assert!(self.contains_q_bucket_idx(hash_bucket_idx));
    //     // TODO: prev bucket id can be out of block scope
    //     if hash_bucket_idx == self.first_q_bucket_idx {
    //         // TODO: optional which is None only on initialization
    //         return self.first_q_bucket_runstart_idx.unwrap();
    //     }
    //     // TODO: no restrictions on hash_bucket_idx but we're expect hash_bucket_idx to belong to
    //     // packed blocks slots
    //     let prev_bucket = hash_bucket_idx.wrapping_sub(1) % self.total_buckets();
    //
    //     (self.run_end(prev_bucket) + 1) % self.total_buckets()
    // }
    //
    pub fn contains_q_bucket_idx(&self, q_bucket_idx: u64) -> bool {
        q_bucket_idx >= self.first_q_bucket_idx && q_bucket_idx < self.first_q_bucket_idx + 64
    }
}

impl Membership for Block {
    fn contains<T: Hash>(&self, _item: T) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_block() {
        let mut f = Filter::new(100, 0.01).unwrap();

        let count_blocks = f.blocks().count();
        assert_eq!(count_blocks, 2);

        let rbits = f.rbits.get() as u8;
        let qbits = f.qbits.get() as u8;

        f.insert_raw_counting(1, (55 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (55 << rbits) | (2u64)).unwrap();
        f.insert_raw_counting(1, (56 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (57 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (57 << rbits) | (2u64)).unwrap();
        f.insert_raw_counting(1, (58 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (2u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (3u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (4u64)).unwrap();
        f.insert_raw_counting(1, (61 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (62 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (62 << rbits) | (3u64)).unwrap();
        f.insert_raw_counting(1, (62 << rbits) | (3u64)).unwrap();
        f.insert_raw_counting(1, (64 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (65 << rbits) | (1u64)).unwrap();

        let raw_block_bytes = f.block_bytes_with_r(0);
        let runstart_idx = f.run_start(0);

        let block = Block {
            first_q_bucket_idx: 0,
            first_q_bucket_runstart_idx: Some(runstart_idx),

            buffer: raw_block_bytes.to_vec(),

            rbits: f.rbits.get() as u8,
            qbits: f.qbits.get() as u8,

            next: None,
        };
        let bucket_idx = 55;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        let bucket_idx = 56;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        let bucket_idx = 57;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        let bucket_idx = 61;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 62;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 63;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
    }

    #[test]
    fn multiple_blocks_with_offset() {
        let mut f = Filter::new(100, 0.01).unwrap();

        let count_blocks = f.blocks().count();
        assert_eq!(count_blocks, 2);

        let rbits = f.rbits.get() as u8;

        f.insert_raw_counting(1, (55 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (55 << rbits) | (2u64)).unwrap();
        f.insert_raw_counting(1, (56 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (57 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (57 << rbits) | (2u64)).unwrap();
        f.insert_raw_counting(1, (58 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (2u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (3u64)).unwrap();
        f.insert_raw_counting(1, (60 << rbits) | (4u64)).unwrap();
        f.insert_raw_counting(1, (61 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (62 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (62 << rbits) | (3u64)).unwrap();
        f.insert_raw_counting(1, (62 << rbits) | (3u64)).unwrap();
        f.insert_raw_counting(1, (64 << rbits) | (1u64)).unwrap();
        f.insert_raw_counting(1, (65 << rbits) | (1u64)).unwrap();

        let raw_block = f.raw_block(0);
        println!("block 0 offset: {}", raw_block.offset);

        let raw_block = f.raw_block(1);
        println!("block 1 offset: {}", raw_block.offset);

        let block_start_bucket_idx = 64;
        let raw_block_bytes = f.block_bytes_with_r(block_start_bucket_idx / 64);
        let runstart_idx = f.run_start(block_start_bucket_idx);
        let block1 = Block {
            first_q_bucket_idx: 64,
            first_q_bucket_runstart_idx: Some(runstart_idx),

            buffer: raw_block_bytes.to_vec(),

            rbits: f.rbits.get() as u8,
            qbits: f.qbits.get() as u8,
            next: None,
        };

        let block_start_bucket_idx = 0;
        let raw_block_bytes = f.block_bytes_with_r(block_start_bucket_idx);
        let runstart_idx = f.run_start(block_start_bucket_idx);
        let block0 = Block {
            first_q_bucket_idx: 0,
            first_q_bucket_runstart_idx: Some(runstart_idx),

            buffer: raw_block_bytes.to_vec(),

            rbits: f.rbits.get() as u8,
            qbits: f.qbits.get() as u8,

            next: Some(Box::new(block1.clone())),
        };

        println!("---------- block 0 -----------");
        let bucket_idx = 55;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block0.run_end(bucket_idx);

        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 56;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block0.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 57;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block0.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 61;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block0.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        let bucket_idx = 62;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block0.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        let bucket_idx = 63;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block0.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        println!("------------------------------");

        println!("---------- block 1 -----------");
        let bucket_idx = 64;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block1.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 65;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block1.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );

        let bucket_idx = 66;
        let f_run_end = f.run_end(bucket_idx);
        let b_run_end = block1.run_end(bucket_idx);
        println!(
            "{} run end: f:{:?}, b:{:?}",
            bucket_idx, f_run_end, b_run_end
        );
        println!("------------------------------");
    }
}
