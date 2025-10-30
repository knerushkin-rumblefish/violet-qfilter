use serde::{Deserialize, Serialize};
use std::hash::{Hash, Hasher};

use crate::{BitExt, CastNonZeroU8, NonZeroU8, StableHasher};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ByteBlock {
    pub offset: usize,
    pub bytes: Vec<u8>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Run {
    pub run_blocks: Option<Vec<ByteBlock>>,
    pub buffer: Vec<u8>,
    pub q_bucket_block: ByteBlock,
    pub q_bucket_idx: u64,
    pub start_idx: Option<u64>,
    pub end_idx: Option<u64>,

    pub qbits: NonZeroU8,
    pub rbits: NonZeroU8,
}

impl Run {
    #[inline]
    pub fn block_byte_size(&self) -> usize {
        1 + 8 + 8 + 64 * self.rbits.usize() / 8
    }

    #[inline]
    fn hash<T: Hash>(&self, item: T) -> u64 {
        let mut hasher = StableHasher::new();
        item.hash(&mut hasher);
        hasher.finish()
    }

    fn is_occupied(&self, q_bucket_idx: u64) -> bool {
        let q_block_bytes = &self.buffer[..1 + 8 + 8];

        let occupieds = u64::from_le_bytes(q_block_bytes[1..][..8].try_into().unwrap());
        occupieds.is_bit_set((q_bucket_idx % 64) as usize)
    }

    #[inline]
    fn calc_qr(&self, hash: u64) -> (u64, u64) {
        let hash_bucket_idx = (hash >> self.rbits.get()) & ((1 << self.qbits.get()) - 1);
        let remainder = hash & ((1 << self.rbits.get()) - 1);
        (hash_bucket_idx, remainder)
    }

    pub fn contains<T: Hash>(&self, item: T) -> bool {
        let hash = self.hash(item);
        let (hash_bucket_idx, hash_remainder) = self.calc_qr(hash);

        let run_blocks = &self.buffer[1 + 8 + 8..];

        let blocks_count = run_blocks.len() / self.block_byte_size();

        let mut is_contains = false;

        if let Some((start_idx, end_idx)) = self.start_idx.zip(self.end_idx) {
            if self.is_occupied(hash_bucket_idx % 64) {
                for block_num in 0..blocks_count {
                    let block_start_bit = block_num * self.block_byte_size();

                    let block_start = if block_num > 0 { block_start_bit } else { 0 };
                    let block_end = block_start + self.block_byte_size();
                    let block_bytes = &run_blocks[block_start..block_end];

                    let reminders = &block_bytes[1 + 8 + 8..][..8 * self.rbits.usize()];

                    let mut run_start_intra_block_idx = 0;
                    let mut run_end_intra_block_idx = 63;
                    // single block
                    if block_num == 0 && blocks_count == 1 {
                        run_start_intra_block_idx = start_idx % 64;
                        run_end_intra_block_idx = end_idx % 64;
                        // first block in block sequence
                    } else if block_num == 0 {
                        run_start_intra_block_idx = start_idx % 64;
                    // last block in block sequence
                    } else if block_num == blocks_count - 1 {
                        run_end_intra_block_idx = end_idx % 64;
                    }
                    for bucket_idx in run_start_intra_block_idx..=run_end_intra_block_idx {
                        let bucket_start_bit = bucket_idx * self.rbits.u64();
                        let bucket_end_bit = bucket_start_bit + self.rbits.u64();

                        let start_u64 = (bucket_start_bit / 64) as usize;
                        let num_rem_parts =
                            1 + (bucket_end_bit > ((start_u64 + 1) * 64) as u64) as usize;

                        let extra_low = bucket_start_bit as usize - start_u64 * 64;
                        let extra_high =
                            ((start_u64 + 1) * 64).saturating_sub(bucket_end_bit as usize);

                        let rem_parts_bytes = &reminders[start_u64 * 8..][..num_rem_parts * 8];
                        let rem_part = u64::from_le_bytes(rem_parts_bytes[..8].try_into().unwrap());
                        let mut remainder = (rem_part << extra_high) >> (extra_high + extra_low);

                        if let Some(rem_part) = rem_parts_bytes.get(8..16) {
                            let remaining_bits = bucket_end_bit - ((start_u64 + 1) * 64) as u64;
                            let rem_part = u64::from_le_bytes(rem_part.try_into().unwrap());
                            remainder |= (rem_part & !(u64::MAX << remaining_bits))
                                << (self.rbits.u64() - remaining_bits);
                        }

                        if hash_remainder == remainder && self.q_bucket_idx == hash_bucket_idx {
                            is_contains = true;
                        }
                    }
                }
            }
        }

        is_contains
    }
}
