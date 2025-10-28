use crate::{BitExt, Block, StableHasher};
use std::hash::{Hash, Hasher};

#[derive(Default, Debug)]
pub struct PackedBlocks {
    pub buffer: Vec<u8>,
    pub rbits: u8,
    pub qbits: u8,
    pub first_block_idx: u64,
    // first run idx in block can be outside of block slots
    pub first_run_idx: Option<u64>,
    // first run block idx is run idx appropriate to first slot of block
    // TODO: should be replaced with offset u8
    pub first_run_block_idx: u64,
    pub blocks_num: u8,
}

impl PackedBlocks {
    // TODO: internal hasher should be replaced on SHA-256 folding to 64 bytes
    #[inline]
    fn hash<T: Hash>(&self, item: T) -> u64 {
        let mut hasher = StableHasher::new();
        item.hash(&mut hasher);
        hasher.finish()
    }

    #[inline]
    pub fn block_byte_size(&self) -> usize {
        1 + 8 + 8 + 64 * self.rbits as usize / 8
    }

    #[inline]
    fn calc_qr(&self, hash: u64) -> (u64, u64) {
        let q_bucket_idx = (hash >> self.rbits) & ((1 << self.qbits) - 1);
        let remainder = hash & ((1 << self.rbits) - 1);
        (q_bucket_idx, remainder)
    }

    pub fn block_bytes(&self, block_num: u64) -> &[u8] {
        // TODO: is it needed on level of pack of blocks. how to behave if last run is spaned
        // across [last, first] blocks? in this case we need to provide multiple DA blocks
        // let block_num = block_num % self.total_blocks();
        let block_offset = block_num as usize * self.block_byte_size();

        &self.buffer[block_offset..][..1 + 8 + 8]
    }

    #[inline]
    fn blocks_count(&self) -> usize {
        self.buffer.len() / self.block_byte_size()
    }

    #[inline]
    fn total_buckets(&self) -> u64 {
        1 << self.qbits
    }

    #[inline]
    pub fn raw_block(&self, block_num: u64) -> Block {
        // no option to have circular block probing
        // let block_num = block_num % self.total_blocks();
        let block_num = self.normilize_block_num(block_num);
        let block_start = block_num as usize * self.block_byte_size();
        println!("PackedBlocks: block_start: {block_start}");
        let block_bytes: &[u8; 1 + 8 + 8] =
            &self.buffer[block_start..][..1 + 8 + 8].try_into().unwrap();
        Block {
            offset: block_bytes[0] as u64,
            occupieds: u64::from_le_bytes(block_bytes[1..1 + 8].try_into().unwrap()),
            runends: u64::from_le_bytes(block_bytes[1 + 8..1 + 8 + 8].try_into().unwrap()),
        }
    }

    /// End idx of the end of the run (inclusive).
    // pub fn run_end(&self, hash_bucket_idx: u64) -> u64 {
    //     let hash_bucket_idx = hash_bucket_idx % self.total_buckets();
    //     let bucket_block_idx = self.normilize_bucket_block_idx(hash_bucket_idx / 64);
    //     let bucket_intrablock_offset = hash_bucket_idx % 64;
    //
    //     let bucket_block = self.block(bucket_block_idx);
    //     let bucket_intrablock_rank = bucket_block.occupieds.popcnt(..=bucket_intrablock_offset);
    //     // No occupied buckets all the way to bucket_intrablock_offset
    //     // which also means hash_bucket_idx isn't occupied
    //     if bucket_intrablock_rank == 0 {
    //         // TODO: what is this?
    //         return if bucket_block.offset <= bucket_intrablock_offset {
    //             // hash_bucket_idx points to an empty bucket unaffected by block offset,
    //             // thus end == start
    //             hash_bucket_idx
    //         } else {
    //             // hash_bucket_idx fall within the section occupied by the offset,
    //             // thus end == last bucket of offset section
    //             (bucket_block_idx * 64 + bucket_block.offset - 1) % self.total_buckets()
    //         };
    //     }
    //
    //     // Must search runends to figure out the end of the run
    //     // println!("run_end bucket block offset: {}", bucket_block.offset);
    //     // Is this block or next one?
    //     let mut runend_block_idx = bucket_block_idx + bucket_block.offset / 64;
    //     // println!("run_end block idx: {runend_block_idx}");
    //     let mut runend_ignore_bits = bucket_block.offset % 64;
    //     let mut runend_block = self.raw_block(runend_block_idx);
    //     // Try to find the runend for the bucket in this block.
    //     // We're looking for the runend_rank'th bit set (0 based)
    //     let mut runend_rank = bucket_intrablock_rank - 1;
    //
    //     // println!("run_end ignore bits: {runend_ignore_bits}");
    //     // println!("run_end rank: {runend_rank}");
    //     let mut runend_block_offset = runend_block
    //         .runends
    //         .select(runend_ignore_bits.., runend_rank);
    //
    //     if let Some(runend_block_offset) = runend_block_offset {
    //         let runend_idx = runend_block_idx * 64 + runend_block_offset;
    //         return runend_idx.max(hash_bucket_idx) % self.total_buckets();
    //     }
    //     // There were not enough runend bits set, keep looking...
    //     loop {
    //         // subtract any runend bits found
    //         runend_rank -= runend_block.runends.popcnt(runend_ignore_bits..);
    //         // println!("run_end rank: {runend_rank}");
    //         // move to the next block
    //         // TODO: should be limited to packed blocks
    //         runend_block_idx += 1;
    //         runend_ignore_bits = 0;
    //         runend_block = self.raw_block(runend_block_idx);
    //         // println!("run_end block: {runend_block_idx}");
    //         runend_block_offset = runend_block
    //             .runends
    //             .select(runend_ignore_bits.., runend_rank);
    //
    //         if let Some(runend_block_offset) = runend_block_offset {
    //             // println!("run_end block offset: {runend_block_offset}");
    //             let runend_idx = runend_block_idx * 64 + runend_block_offset;
    //             // println!("run_end runend_idx: {runend_idx}");
    //             // println!("run_end total buckets: {:}", self.total_buckets());
    //             //
    //             // println!(
    //             //     "run_end runend_idx.max(hash_bucket_idx): {:}",
    //             //     runend_idx.max(hash_bucket_idx)
    //             // );
    //             // println!(
    //             //     "run_end runend_idx.max(hash_bucket_idx) % self.total_buckets(): {:}",
    //             //     runend_idx.max(hash_bucket_idx) % self.total_buckets(),
    //             // );
    //             return runend_idx.max(hash_bucket_idx) % self.total_buckets();
    //         }
    //     }
    // }

    // #[cold]
    // #[inline(never)]
    // fn calc_offset(&self, block_num: u64) -> u64 {
    //     // TODO: should buffer on packed blocks level contains padding to 482 bytes?
    //     // the block offset can be calculated as the difference between its position and runstart.
    //     let block_start = block_num * 64;
    //     let mut run_start = self.run_start(block_start);
    //     // TODO: how to handle this rotating property of linear probing
    //     // in regular case run_start is >= block_start but if it's rotate over buffer end its
    //     // circulate in beginning of buffer, in case of packing it by chunks it's require to
    //     // provide another chunk
    //     // if run_start < block_start {
    //     //     run_start += self.total_buckets().get();
    //     // }
    //     run_start - block_start
    // }

    // #[inline]
    // pub fn block(&self, block_num: u64) -> Block {
    //     let block_bytes: &[u8; 1 + 8 + 8] = &self.block_bytes(block_num).try_into().unwrap();
    //     let offset = {
    //         if block_bytes[0] < u8::MAX {
    //             block_bytes[0] as u64
    //         } else {
    //             self.calc_offset(block_num)
    //         }
    //     };
    //     Block {
    //         offset,
    //         occupieds: u64::from_le_bytes(block_bytes[1..1 + 8].try_into().unwrap()),
    //         runends: u64::from_le_bytes(block_bytes[1 + 8..1 + 8 + 8].try_into().unwrap()),
    //     }
    // }

    pub fn validate_bucket(&self, bucket_block_idx: u64) -> bool {
        bucket_block_idx >= self.first_block_idx
            && bucket_block_idx
                < self.first_block_idx
                    + ((self.blocks_num as u64) * (self.block_byte_size() as u64))
    }

    pub fn normilize_bucket_block_idx(&self, bucket_block_idx: u64) -> u64 {
        assert!(self.validate_bucket(bucket_block_idx));

        println!(
            "PackedBlocks: normilize_bucket_block_idx: bucket_block_idx: {}",
            bucket_block_idx - self.first_block_idx
        );
        bucket_block_idx - self.first_block_idx
    }

    pub fn normilize_block_num(&self, block_num: u64) -> u64 {
        println!("PackedBlocks: normilize_block_num: block_num: {block_num}");
        (self.normilize_bucket_block_idx(block_num * 64) + 1) / 64
    }

    // pub fn contains<T: Hash>(&self, item: T) -> bool {
    //     let hash = self.hash(item);
    //     let (q_bucket_idx, remainder) = self.calc_qr(hash);
    //
    //     let bucket_block_idx = self.normilize_bucket_block_idx(q_bucket_idx / 64);
    //     let bucket_intrablock_offset = q_bucket_idx % 64;
    //
    //     let bucket_block = self.block(bucket_block_idx);
    // }
}
