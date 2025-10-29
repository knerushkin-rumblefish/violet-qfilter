use crate::{BitExt, Block, StableHasher};
use std::hash::{Hash, Hasher};

#[derive(Default, Debug)]
pub struct PackedBlocks {
    pub buffer: Vec<u8>,
    pub rbits: u8,
    pub qbits: u8,
    pub block_offset: u64,
    pub blocks_num: u8,

    pub first_runstart_idx: u64,
    //          block0 [
    //                   (run02 (s02N)],
    //                   [run01 (s011, s012, .., s011N)],
    //                 ],
    //        ________________________________________________
    // packed | block1 [                                     |
    //        |           [run11 (s111, s112, .., s11N)]     |
    //        |        ],                                    |
    //        | block2 [                                     |
    //        |          [run21 (s211, s212, .., s21N)],     |
    //    >>> |          [run22 (s221, s222, .., s22(N - 1)) | >>> shifted
    //        |        ],                                    |
    //        ________________________________________________
    //          block3 [
    //                   (run22 (s22N)],
    //                   [run31 (s311, s312, .., s311N)],
    //                 ],
    pub shifted: bool,
}

impl PackedBlocks {
    // in case PackedBlocks not contain all runs slots in buffer it should require next PackedBlock
    // extension do not extened range of runend's beloning to PackedBlocks, its just extend memory
    // to run start belonging to next block so we can call contain() on remainders belonding to block slots
    pub fn extend(&mut self, extension: &mut PackedBlocks) {
        self.buffer.append(&mut extension.buffer);
        // should be just calculated from buffer size
        self.blocks_num += extension.blocks_num;
        // TODO: extension should cover all block slots remainders
        // so extension should be validated that any remainder belonging to block is covered by
        // current buffer and not shifted to next one
        // in case of membership proof it's can be enough to have croped memory, but in
        // case of non-membership we should traverse all valued corresponding to run
    }

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

    fn block_bytes(&self, block_num: u64) -> &[u8] {
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
        //
        let block_num = self.normilize_block_num(block_num);
        let block_start = block_num as usize * self.block_byte_size();
        let block_bytes: &[u8; 1 + 8 + 8] =
            &self.buffer[block_start..][..1 + 8 + 8].try_into().unwrap();
        Block {
            offset: block_bytes[0] as u64,
            occupieds: u64::from_le_bytes(block_bytes[1..1 + 8].try_into().unwrap()),
            runends: u64::from_le_bytes(block_bytes[1 + 8..1 + 8 + 8].try_into().unwrap()),
        }
    }

    // End idx of the end of the run (inclusive).
    pub fn run_end(&self, hash_bucket_idx: u64) -> u64 {
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();
        let bucket_block_idx = hash_bucket_idx / 64;
        let bucket_intrablock_offset = hash_bucket_idx % 64;

        let bucket_block = self.block(bucket_block_idx);
        let bucket_intrablock_rank = bucket_block.occupieds.popcnt(..=bucket_intrablock_offset);
        // No occupied buckets all the way to bucket_intrablock_offset
        // which also means hash_bucket_idx isn't occupied
        if bucket_intrablock_rank == 0 {
            // TODO: WTF???
            return if bucket_block.offset <= bucket_intrablock_offset {
                // hash_bucket_idx points to an empty bucket unaffected by block offset,
                // thus end == start
                hash_bucket_idx
            } else {
                // hash_bucket_idx fall within the section occupied by the offset,
                // thus end == last bucket of offset section
                // TODO: WTF???
                (bucket_block_idx * 64 + bucket_block.offset - 1) % self.total_buckets()
            };
        }

        // Must search runends to figure out the end of the run
        // println!("run_end bucket block offset: {}", bucket_block.offset);
        // Is this block or next one?
        let mut runend_block_idx = bucket_block_idx + bucket_block.offset / 64;
        // println!("run_end block idx: {runend_block_idx}");
        let mut runend_ignore_bits = bucket_block.offset % 64;
        let mut runend_block = self.raw_block(runend_block_idx);
        // Try to find the runend for the bucket in this block.
        // We're looking for the runend_rank'th bit set (0 based)
        let mut runend_rank = bucket_intrablock_rank - 1;

        // println!("run_end ignore bits: {runend_ignore_bits}");
        // println!("run_end rank: {runend_rank}");
        let mut runend_block_offset = runend_block
            .runends
            .select(runend_ignore_bits.., runend_rank);

        if let Some(runend_block_offset) = runend_block_offset {
            let runend_idx = runend_block_idx * 64 + runend_block_offset;
            return runend_idx.max(hash_bucket_idx) % self.total_buckets();
        }
        // There were not enough runend bits set, keep looking...
        loop {
            // subtract any runend bits found
            runend_rank -= runend_block.runends.popcnt(runend_ignore_bits..);
            // println!("run_end rank: {runend_rank}");
            // move to the next block
            // TODO: should be limited to packed blocks
            runend_block_idx += 1;
            runend_ignore_bits = 0;
            runend_block = self.raw_block(runend_block_idx);
            // println!("run_end block: {runend_block_idx}");
            runend_block_offset = runend_block
                .runends
                .select(runend_ignore_bits.., runend_rank);

            if let Some(runend_block_offset) = runend_block_offset {
                // println!("run_end block offset: {runend_block_offset}");
                let runend_idx = runend_block_idx * 64 + runend_block_offset;
                // println!("run_end runend_idx: {runend_idx}");
                // println!("run_end total buckets: {:}", self.total_buckets());
                //
                // println!(
                //     "run_end runend_idx.max(hash_bucket_idx): {:}",
                //     runend_idx.max(hash_bucket_idx)
                // );
                // println!(
                //     "run_end runend_idx.max(hash_bucket_idx) % self.total_buckets(): {:}",
                //     runend_idx.max(hash_bucket_idx) % self.total_buckets(),
                // );
                return runend_idx.max(hash_bucket_idx) % self.total_buckets();
            }
        }
    }

    #[inline]
    fn run_start(&self, hash_bucket_idx: u64) -> u64 {
        // TODO: prev bucket id can be out of packed blocs scope
        let first_block_first_bucket_idx = self.block_offset * 64;
        if hash_bucket_idx == first_block_first_bucket_idx {
            return self.first_runstart_idx;
        }
        // TODO: no restrictions on hash_bucket_idx but we're expect hash_bucket_idx to belong to
        // packed blocks slots
        let prev_bucket = hash_bucket_idx.wrapping_sub(1) % self.total_buckets();

        (self.run_end(prev_bucket) + 1) % self.total_buckets()
    }

    #[cold]
    #[inline(never)]
    fn calc_offset(&self, block_num: u64) -> u64 {
        // TODO: should buffer on packed blocks level contains padding to 482 bytes?
        // the block offset can be calculated as the difference between its position and runstart.
        let block_start = block_num * 64;
        let run_start = self.run_start(block_start);
        // TODO: how to handle this rotating property of linear probing
        // in regular case run_start is >= block_start but if it's rotate over buffer end its
        // circulate in beginning of buffer, in case of packing it by chunks it's require to
        // provide another chunk
        // if run_start < block_start {
        //     run_start += self.total_buckets().get();
        // }
        run_start - block_start
    }

    #[inline]
    pub fn block(&self, block_num: u64) -> Block {
        let block_num = self.normilize_block_num(block_num);
        let block_bytes: &[u8; 1 + 8 + 8] = &self.block_bytes(block_num).try_into().unwrap();
        let offset = {
            if block_bytes[0] < u8::MAX {
                block_bytes[0] as u64
            } else {
                self.calc_offset(block_num)
            }
        };
        Block {
            offset,
            occupieds: u64::from_le_bytes(block_bytes[1..1 + 8].try_into().unwrap()),
            runends: u64::from_le_bytes(block_bytes[1 + 8..1 + 8 + 8].try_into().unwrap()),
        }
    }

    pub fn validate_bucket_idx(&self, bucket_idx: u64) -> bool {
        // block 1 [idx1, idx2, ..., idxN], ..., block N [idx1, idx2, ..., idxN]
        //          ^first_block_first_bucket_idx                          ^last_block_last_bucket_idx
        let first_block_first_bucket_idx = self.block_offset * 64;
        let last_block_last_bucket_idx =
            first_block_first_bucket_idx + ((self.blocks_num as u64) * 64) - 1;
        bucket_idx >= first_block_first_bucket_idx && bucket_idx < last_block_last_bucket_idx
    }

    pub fn normilize_bucket_idx(&self, bucket_idx: u64) -> u64 {
        assert!(self.validate_bucket_idx(bucket_idx));
        let first_block_first_bucket_idx = self.block_offset * 64;
        bucket_idx - first_block_first_bucket_idx
    }

    pub fn normilize_block_num(&self, block_num: u64) -> u64 {
        let block_bucket_idx = block_num * 64;
        // [block1][block2][block3][block4][block5]...
        //         [^^^^^^^packed^^^^^^^^^]
        //         [block0][block1][block2]
        (self.normilize_bucket_idx(block_bucket_idx) + 1) / 64
    }

    pub fn runs_span(&self) -> (u64, u64) {
        //TODO: should return (first, last) run idx contained in packed blocks and belondign to
        //slot of packed blocks
        (
            (self.block_offset * 64),
            (self.block_offset * 64) + (self.blocks_num as u64 * 64) - 1,
        )
    }

    #[inline(always)]
    fn is_occupied(&self, hash_bucket_idx: u64) -> bool {
        println!(
            "packed blocks: is_occupied hash_bukcet_idx {}",
            hash_bucket_idx
        );
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();
        let hash_bucket_block_idx = hash_bucket_idx / 64;
        let hash_bucket_block_idx = self.normilize_block_num(hash_bucket_block_idx);
        println!(
            "packed blocks: is_occupied block_num {}",
            hash_bucket_block_idx
        );
        let block_start = hash_bucket_block_idx as usize * self.block_byte_size();
        println!("packed blocks: is_occupied block_start {}", block_start);
        let occupieds = u64::from_le_bytes(self.buffer[block_start + 1..][..8].try_into().unwrap());
        occupieds.is_bit_set((hash_bucket_idx % 64) as usize)
    }

    #[inline(always)]
    fn is_runend(&self, hash_bucket_idx: u64) -> bool {
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();
        let hash_bucket_block_idx = hash_bucket_idx / 64;
        let hash_bucket_block_idx = self.normilize_block_num(hash_bucket_block_idx);
        let block_start = hash_bucket_block_idx as usize * self.block_byte_size();

        let runends =
            u64::from_le_bytes(self.buffer[block_start + 1 + 8..][..8].try_into().unwrap());
        runends.is_bit_set((hash_bucket_idx % 64) as usize)
    }

    #[inline(always)]
    pub fn get_remainder(&self, hash_bucket_idx: u64) -> u64 {
        println!("packed blocks: get_remainder");
        debug_assert!(self.rbits > 0 && self.rbits < 64);
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();

        let hash_bucket_block_idx = hash_bucket_idx / 64;
        let hash_bucket_block_idx = self.normilize_block_num(hash_bucket_block_idx);

        println!(
            "packed blocks: get_remainder hash_bucket_block_idx: {}",
            hash_bucket_block_idx
        );
        let block_start = hash_bucket_block_idx as usize * self.block_byte_size();

        println!("packed blocks: get_remainder block_start: {}", block_start);

        let remainders_start = block_start + 1 + 8 + 8;
        let start_bit_idx = self.rbits as usize * (hash_bucket_idx % 64) as usize;
        let end_bit_idx = start_bit_idx + self.rbits as usize;
        let start_u64 = start_bit_idx / 64;
        let num_rem_parts = 1 + (end_bit_idx > (start_u64 + 1) * 64) as usize;
        let rem_parts_bytes = &self.buffer[remainders_start + start_u64 * 8..][..num_rem_parts * 8];
        let extra_low = start_bit_idx - start_u64 * 64;
        let extra_high = ((start_u64 + 1) * 64).saturating_sub(end_bit_idx);
        let rem_part = u64::from_le_bytes(rem_parts_bytes[..8].try_into().unwrap());
        // zero high bits & truncate low bits
        let mut remainder = (rem_part << extra_high) >> (extra_high + extra_low);
        if let Some(rem_part) = rem_parts_bytes.get(8..16) {
            let remaining_bits = end_bit_idx - (start_u64 + 1) * 64;
            let rem_part = u64::from_le_bytes(rem_part.try_into().unwrap());
            remainder |= (rem_part & !(u64::MAX << remaining_bits))
                << (self.rbits as usize - remaining_bits);
        }
        debug_assert!(remainder.leading_zeros() >= 64 - self.rbits as u32);
        remainder
    }

    pub fn contains<T: Hash>(&self, item: T) -> bool {
        let hash = self.hash(item);
        let (hash_bucket_idx, hash_remainder) = self.calc_qr(hash);
        println!(
            "packed blocks contains: hash_bucket_idx {}",
            hash_bucket_idx
        );

        if !self.validate_bucket_idx(hash_bucket_idx) {
            return false;
        }

        if !self.is_occupied(hash_bucket_idx) {
            return false;
        }

        let mut runstart_idx = self.run_start(hash_bucket_idx);

        loop {
            println!("packed blocks contains: runstart_idx: {}", runstart_idx);
            let reminder = self.get_remainder(runstart_idx);
            println!("packed blocks contains: reminder: {}", reminder);
            println!("packed blocks contains: hash_reminder: {}", hash_remainder);
            if hash_remainder == reminder {
                return true;
            }
            if self.is_runend(runstart_idx) {
                return false;
            }
            runstart_idx += 1;
        }
    }
}
