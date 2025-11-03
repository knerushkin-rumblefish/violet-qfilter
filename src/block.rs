use crate::{BitExt, Block as BytesBlock, Filter, Membership, StableHasher};

use std::hash::{Hash, Hasher};

const BLOCK_OFFSET_SIZE: usize = 1;
const BLOCK_OCCUPIEDS_SIZE: usize = 8;
const BLOCK_RUNENDS_SIZE: usize = 8;

const BLOCK_REMINDERS_SLOTS_SIZE: usize = 64;

#[derive(Default, Debug, Clone)]
pub struct Block {
    pub first_q_bucket_idx: u64,
    pub first_q_bucket_runstart_offset: Option<u8>,

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

    #[inline]
    fn hash<T: Hash>(&self, item: T) -> u64 {
        let mut hasher = StableHasher::new();
        item.hash(&mut hasher);
        hasher.finish()
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

        if !self.contains_q_bucket_idx(hash_bucket_idx) {
            return None;
        }

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

        let mut runend_block_idx = hash_bucket_block_idx + current_block.offset / 64;
        // run_end shift should be bounded with single block
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
            } else {
                return None;
            }
        }
    }

    #[inline]
    fn run_start(&self, hash_bucket_idx: u64) -> Option<u64> {
        if !self.contains_q_bucket_idx(hash_bucket_idx) {
            return None;
        }

        // TODO: prev bucket id can be out of block scope
        if hash_bucket_idx == self.first_q_bucket_idx {
            // TODO: optional which is None only on initialization
            return self
                .first_q_bucket_runstart_offset
                .map(|runstart_offset| runstart_offset as u64 + self.first_q_bucket_idx);
        }
        // TODO: no restrictions on hash_bucket_idx but we're expect hash_bucket_idx to belong to
        // packed blocks slots
        let prev_bucket = hash_bucket_idx.wrapping_sub(1) % self.total_buckets();

        if let Some(run_end) = self.run_end(prev_bucket) {
            return Some((run_end + 1) % self.total_buckets());
        }

        None
    }

    #[inline(always)]
    pub fn get_remainder(&self, hash_bucket_idx: u64) -> Option<u64> {
        debug_assert!(self.rbits > 0 && self.rbits < 64);
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();
        if !self.contains_q_bucket_idx(hash_bucket_idx) {
            if let Some(next_block) = &self.next {
                return next_block.get_remainder(hash_bucket_idx);
            }
            return None;
        }

        let start_bit_idx = self.rbits as usize * (hash_bucket_idx % 64) as usize;
        let end_bit_idx = start_bit_idx + self.rbits as usize;
        let start_u64 = start_bit_idx / 64;
        let num_rem_parts = 1 + (end_bit_idx > (start_u64 + 1) * 64) as usize;
        let rem_parts_bytes =
            &self.buffer[self.block_header_size() + start_u64 * 8..][..num_rem_parts * 8];
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
        Some(remainder)
    }

    #[inline(always)]
    fn is_occupied(&self, hash_bucket_idx: u64) -> Option<bool> {
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();

        if !self.contains_q_bucket_idx(hash_bucket_idx) {
            return None;
        }

        let occupieds = u64::from_le_bytes(
            self.buffer[BLOCK_OFFSET_SIZE..][..BLOCK_OCCUPIEDS_SIZE]
                .try_into()
                .unwrap(),
        );
        Some(occupieds.is_bit_set((hash_bucket_idx % 64) as usize))
    }

    #[inline(always)]
    fn is_runend(&self, hash_bucket_idx: u64) -> Option<bool> {
        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();

        if !self.contains_q_bucket_idx(hash_bucket_idx) {
            if let Some(next_block) = &self.next {
                return next_block.is_runend(hash_bucket_idx);
            }
            return None;
        }

        let runends = u64::from_le_bytes(
            self.buffer[BLOCK_OFFSET_SIZE + BLOCK_OCCUPIEDS_SIZE..][..BLOCK_RUNENDS_SIZE]
                .try_into()
                .unwrap(),
        );
        Some(runends.is_bit_set((hash_bucket_idx % 64) as usize))
    }

    pub fn contains_fingerprint(&self, hash: u64) -> Option<bool> {
        let (hash_bucket_idx, hash_remainder) = self.calc_qr(hash);

        let hash_bucket_idx = hash_bucket_idx % self.total_buckets();

        if !self.contains_q_bucket_idx(hash_bucket_idx) {
            return None;
        }

        if !self.is_occupied(hash_bucket_idx)? {
            return None;
        }

        let mut runstart_idx = self.run_start(hash_bucket_idx)?;

        loop {
            let remainder = self.get_remainder(runstart_idx)?;
            if hash_remainder == remainder {
                return Some(true);
            }
            if self.is_runend(runstart_idx)? {
                return Some(false);
            }
            runstart_idx += 1;
        }
    }

    pub fn contains_q_bucket_idx(&self, q_bucket_idx: u64) -> bool {
        q_bucket_idx >= self.first_q_bucket_idx && q_bucket_idx < self.first_q_bucket_idx + 64
    }

    pub fn extract_block_by_fingerprint(filter: &Filter, hash: u64) -> Block {
        let (hash_bucket_idx, _hash_remainder) = filter.calc_qr(hash);

        let hash_bucket_idx = hash_bucket_idx % filter.total_buckets();
        let hash_bucket_block_num = hash_bucket_idx / 64;
        let hash_bucket_block_idx = (hash_bucket_block_num * 64) % filter.total_buckets();

        let raw_block_bytes = filter.block_bytes_with_r(hash_bucket_block_num);

        let hash_bucket_runstart_idx = filter.run_start(hash_bucket_idx);
        let hash_bucket_runend_idx = filter.run_end(hash_bucket_idx);

        // TODO: case when first bucket is already not in the current block
        let first_bucket_runstart_idx = filter.run_start(hash_bucket_block_idx);

        let mut next_block: Option<Box<Block>> = None;
        let next_block_bucket_idx = (hash_bucket_block_idx + 64) % filter.total_buckets();

        // TODO: run shifted more than to the next block
        if hash_bucket_runend_idx >= next_block_bucket_idx
            || hash_bucket_runstart_idx >= next_block_bucket_idx
            || first_bucket_runstart_idx > next_block_bucket_idx
        {
            let next_block_first_bucket_runstart_idx = filter.run_start(next_block_bucket_idx);
            let raw_block_bytes = filter.block_bytes_with_r(next_block_bucket_idx / 64);
            next_block = Some(Box::new(Block {
                first_q_bucket_idx: next_block_bucket_idx,
                first_q_bucket_runstart_offset: Some(
                    (next_block_first_bucket_runstart_idx - next_block_bucket_idx)
                        .try_into()
                        .unwrap(),
                ),

                buffer: raw_block_bytes.to_vec(),

                rbits: filter.rbits.get() as u8,
                qbits: filter.qbits.get() as u8,

                // TODO: should link all following block's on which current run is spanned
                next: None,
            }))
        }

        Block {
            first_q_bucket_idx: hash_bucket_block_idx,
            first_q_bucket_runstart_offset: Some(
                (first_bucket_runstart_idx - hash_bucket_block_idx)
                    .try_into()
                    .unwrap(),
            ),

            buffer: raw_block_bytes.to_vec(),

            rbits: filter.rbits.get() as u8,
            qbits: filter.qbits.get() as u8,

            // TODO: should link all following block's on which current run is spanned
            next: next_block,
        }
    }

    pub fn runs(&self) -> RunIter<'_> {
        RunIter::new(self)
    }
}

pub struct RunIter<'a> {
    block: &'a Block,
    q_bucket_idx: u64,
    r_bucket_idx: u64,
    remaining: u64,
}

impl<'a> RunIter<'a> {
    pub fn new(block: &'a Block) -> Self {
        let mut iter = RunIter {
            block,
            q_bucket_idx: 0,
            r_bucket_idx: 0,
            remaining: 64,
        };

        while let Some(false) = block.is_occupied(iter.q_bucket_idx) {
            iter.q_bucket_idx += 1;
        }
        if let Some(run_start) = block.run_start(iter.q_bucket_idx) {
            iter.r_bucket_idx = run_start;
        }

        iter
    }
}

impl<'a> Iterator for RunIter<'a> {
    type Item = (u64, u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(r) = self.remaining.checked_sub(1) {
            self.remaining = r;
        } else {
            return None;
        }

        let run_start_idx = self.block.run_start(self.q_bucket_idx)?;
        let run_end_idx = self.block.run_end(self.q_bucket_idx)?;

        let run = (self.q_bucket_idx, run_start_idx, run_end_idx);

        if self.block.is_runend(self.r_bucket_idx)? {
            self.q_bucket_idx += 1;
            while let Some(false) = self.block.is_occupied(self.q_bucket_idx) {
                self.q_bucket_idx += 1;
            }
            self.r_bucket_idx = (self.r_bucket_idx + 1).max(self.q_bucket_idx);
        } else {
            self.r_bucket_idx += 1;
        }

        Some(run)
    }
}

impl Membership for Block {
    fn contains<T: Hash>(&self, item: T) -> Option<bool> {
        self.contains_fingerprint(self.hash(item))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Filter;
    use rand::{seq::IndexedRandom, Rng};

    fn generate_filter(n: u64, slot_layout: &[u64]) -> (Filter, Vec<u64>, Vec<u64>) {
        let mut f = Filter::new(n, 0.01).unwrap();

        let count_blocks = f.blocks().count();
        assert_eq!(count_blocks, 2);

        let rbits = f.rbits.get() as u8;
        let qbits = f.qbits.get() as u8;

        let mut remainders: Vec<u64> = vec![];
        let mut fingerprints: Vec<u64> = vec![];

        for &slot in slot_layout.iter() {
            let (mut quotient, mut remainder) = generate_item(slot, qbits, rbits);
            if remainders.contains(&remainder) {
                (quotient, remainder) = generate_item(slot, qbits, rbits);
            }
            remainders.push(remainder);
            fingerprints.push(quotient | remainder);
            // TODO: handle collision. currently expecting remainder to be unique
            f.insert_raw_counting(1, quotient | remainder).unwrap();
        }

        (f, remainders, fingerprints)
    }

    fn generate_item(q_bucket_idx: u64, _qbits: u8, rbits: u8) -> (u64, u64) {
        let mut rng = rand::rng();
        let quotient = q_bucket_idx << rbits;
        let remainder = rng.random::<u64>() & ((1 << rbits) - 1);

        (quotient, remainder)
    }

    fn extract_block(f: &Filter, n: u64, next: Option<Box<Block>>) -> Block {
        let raw_block_bytes = f.block_bytes_with_r(n);

        let first_q_bucket_idx = n * 64;
        let first_q_bucket_runstart_idx = f.run_start(first_q_bucket_idx);

        let block = Block {
            first_q_bucket_idx,
            first_q_bucket_runstart_offset: Some(
                (first_q_bucket_runstart_idx - first_q_bucket_idx)
                    .try_into()
                    .unwrap(),
            ),

            buffer: raw_block_bytes.to_vec(),

            rbits: f.rbits.get() as u8,
            qbits: f.qbits.get() as u8,

            next,
        };

        return block;
    }

    #[test]
    fn single_block_no_next_run_ends() {
        let buckets_layout = vec![55, 55, 56, 57, 57, 58, 59, 60, 60, 60, 61, 62, 62, 64, 65];

        let expected_block_run_ends = vec![56, 56, 57, 59, 59, 60, 61];

        let (f, _remainders, _) = generate_filter(100, &buckets_layout);

        let block = extract_block(&f, 0, None);

        let mut actual_block_run_ends = vec![];
        for bucket_idx in buckets_layout {
            let f_run_end = f.run_end(bucket_idx);
            if let Some(b_run_end) = block.run_end(bucket_idx) {
                actual_block_run_ends.push(b_run_end);
                assert_eq!(f_run_end, b_run_end);
            }
        }

        assert_eq!(expected_block_run_ends, actual_block_run_ends);
    }

    #[test]
    fn multiple_blocks_with_next_run_ends() {
        let buckets_layout = vec![
            55, 55, 56, 57, 57, 58, 59, 60, 60, 60, 61, 62, 62, 63, 64, 64, 65,
        ];

        let (f, _remainders, _) = generate_filter(100, &buckets_layout);
        let expected_block_run_ends = vec![56, 56, 57, 59, 59, 60, 61, 64, 64, 64, 65, 67, 67, 68];

        let block_next = extract_block(&f, 1, None);
        let next = Some(Box::new(block_next));

        let block = extract_block(&f, 0, next);

        let mut actual_block_run_ends = vec![];

        let block_bucket_range = block.first_q_bucket_idx..block.first_q_bucket_idx + 64;
        for bucket_idx in buckets_layout {
            if block_bucket_range.contains(&bucket_idx) {
                let f_run_end = f.run_end(bucket_idx);
                if let Some(b_run_end) = block.run_end(bucket_idx) {
                    actual_block_run_ends.push(b_run_end);
                    assert_eq!(f_run_end, b_run_end);
                }
            }
        }

        assert_eq!(expected_block_run_ends, actual_block_run_ends);
    }

    #[test]
    fn single_block_no_next_run_starts() {
        let buckets_layout = vec![55, 55, 56, 57, 57, 58, 59, 60, 60, 60, 61, 62, 62, 64, 65];

        let expected_block_run_starts = vec![55, 55, 57, 58, 58, 60, 61, 62, 62, 62];

        let (f, _remainders, _) = generate_filter(100, &buckets_layout);

        let next = None;

        let block = extract_block(&f, 0, next);

        let mut actual_block_run_starts = vec![];

        for bucket_idx in buckets_layout {
            let f_run_start = f.run_start(bucket_idx);
            if let Some(b_run_start) = block.run_start(bucket_idx) {
                actual_block_run_starts.push(b_run_start);
                assert_eq!(f_run_start, b_run_start);
            }
        }

        assert_eq!(actual_block_run_starts, expected_block_run_starts);
    }

    #[test]
    fn multiple_blocks_with_next_run_starts() {
        let buckets_layout = vec![55, 55, 56, 57, 57, 58, 59, 60, 60, 60, 61, 62, 62, 64, 65];

        let expected_block_run_starts = vec![55, 55, 57, 58, 58, 60, 61, 62, 62, 62, 65, 66, 66];

        let (f, _remainders, _) = generate_filter(100, &buckets_layout);

        let block_next = extract_block(&f, 1, None);
        let next = Some(Box::new(block_next));

        let block = extract_block(&f, 0, next);
        let mut actual_block_run_starts = vec![];

        for bucket_idx in buckets_layout {
            let f_run_start = f.run_start(bucket_idx);
            if let Some(b_run_start) = block.run_start(bucket_idx) {
                actual_block_run_starts.push(b_run_start);
                assert_eq!(f_run_start, b_run_start);
            }
        }

        assert_eq!(actual_block_run_starts, expected_block_run_starts);
    }

    #[test]
    fn single_block_no_next_no_collision_remainders() {
        let buckets_layout: Vec<u64> = (55..66).collect();

        let (f, remainders, _) = generate_filter(100, &buckets_layout);

        let expected_block_remainders: &[u64] = &remainders[..buckets_layout
            .clone()
            .into_iter()
            .filter(|&bucket_idx| bucket_idx < 64)
            .count()];

        let next = None;

        let block = extract_block(&f, 0, next);

        let mut actual_block_remainders = vec![];

        for bucket_idx in buckets_layout {
            let f_remainder = f.get_remainder(bucket_idx);
            if let Some(b_remainder) = block.get_remainder(bucket_idx) {
                actual_block_remainders.push(b_remainder);
                assert_eq!(f_remainder, b_remainder);
            }
        }

        assert_eq!(actual_block_remainders, expected_block_remainders);
    }

    #[test]
    fn multiple_block_with_next_no_collision_remainders() {
        let buckets_layout: Vec<u64> = (55..66).collect();

        let (f, remainders, _) = generate_filter(100, &buckets_layout);

        let expected_block_remainders: &[u64] = &remainders[..buckets_layout
            .clone()
            .into_iter()
            .filter(|&bucket_idx| bucket_idx < 64)
            .count()];

        let block_next = extract_block(&f, 1, None);
        let next = Some(Box::new(block_next));

        let block = extract_block(&f, 0, next);

        let mut actual_block_remainders = vec![];

        for bucket_idx in buckets_layout {
            let f_run_end = f.run_end(bucket_idx);
            let f_remainder = f.get_remainder(f_run_end);
            if let Some(b_run_end) = block.run_end(bucket_idx) {
                if let Some(b_remainder) = block.get_remainder(b_run_end) {
                    actual_block_remainders.push(b_remainder);
                    assert_eq!(f_remainder, b_remainder);
                }
            }
        }
        assert_eq!(actual_block_remainders, expected_block_remainders);
    }

    #[test]
    fn multiple_block_with_next_single_collision_remainders() {
        let mut buckets_layout: Vec<u64> = vec![0; 10];
        buckets_layout[0] = 57;
        buckets_layout[1] = 58;
        buckets_layout[2] = 58;
        buckets_layout[3] = 59;
        buckets_layout[4] = 60;
        buckets_layout[5] = 61;
        buckets_layout[6] = 62;
        buckets_layout[7] = 63;
        buckets_layout[8] = 64;
        buckets_layout[9] = 65;

        let (f, remainders, _) = generate_filter(100, &buckets_layout);

        let mut collision_remainders = [remainders[1], remainders[2]];
        collision_remainders.sort();

        let expected_reminders_size = buckets_layout
            .iter()
            .filter(|&bucket_idx| *bucket_idx < 64)
            .count();

        let mut expected_block_remainders: Vec<u64> = vec![0; 10];
        expected_block_remainders[0] = remainders[0];
        expected_block_remainders[1] = collision_remainders[1];
        expected_block_remainders[2] = collision_remainders[1];
        expected_block_remainders[3] = remainders[3];
        expected_block_remainders[4] = remainders[4];
        expected_block_remainders[5] = remainders[5];
        expected_block_remainders[6] = remainders[6];
        expected_block_remainders[7] = remainders[7];
        expected_block_remainders[8] = remainders[8];
        expected_block_remainders[9] = remainders[9];
        let expected_block_remainders = &expected_block_remainders[..expected_reminders_size];

        let block_next = extract_block(&f, 1, None);
        let next = Some(Box::new(block_next));

        let block = extract_block(&f, 0, next);

        let mut actual_block_remainders = vec![];

        for bucket_idx in buckets_layout {
            let f_run_end = f.run_end(bucket_idx);
            let f_remainder = f.get_remainder(f_run_end);
            if let Some(b_run_end) = block.run_end(bucket_idx) {
                if let Some(b_remainder) = block.get_remainder(b_run_end) {
                    actual_block_remainders.push(b_remainder);
                    assert_eq!(f_remainder, b_remainder);
                }
            }
        }

        assert_eq!(actual_block_remainders, expected_block_remainders);
    }

    #[test]
    fn single_block_no_next_single_collision_remainders() {
        let mut buckets_layout: Vec<u64> = vec![0; 10];
        buckets_layout[0] = 57;
        buckets_layout[1] = 58;
        buckets_layout[2] = 58;
        buckets_layout[3] = 59;
        buckets_layout[4] = 60;
        buckets_layout[5] = 61;
        buckets_layout[6] = 62;
        buckets_layout[7] = 63;
        buckets_layout[8] = 64;
        buckets_layout[9] = 65;

        let expected_run_end = [57, 59, 59, 60, 61, 62, 63, 64, 65, 66];
        let (f, remainders, _) = generate_filter(100, &buckets_layout);

        let mut collision_remainders = [remainders[1], remainders[2]];
        collision_remainders.sort();

        let expected_reminders_size = expected_run_end
            .iter()
            .filter(|&bucket_idx| *bucket_idx < 64)
            .count();

        let mut expected_block_remainders: Vec<u64> = vec![0; 10];
        expected_block_remainders[0] = remainders[0];
        expected_block_remainders[1] = collision_remainders[1];
        expected_block_remainders[2] = collision_remainders[1];
        expected_block_remainders[3] = remainders[3];
        expected_block_remainders[4] = remainders[4];
        expected_block_remainders[5] = remainders[5];
        expected_block_remainders[6] = remainders[6];
        expected_block_remainders[7] = remainders[7];
        expected_block_remainders[8] = remainders[8];
        expected_block_remainders[9] = remainders[9];
        let expected_block_remainders = &expected_block_remainders[..expected_reminders_size];

        let next = None;

        let block = extract_block(&f, 0, next);

        let mut actual_block_remainders = vec![];

        for bucket_idx in buckets_layout {
            let f_run_end = f.run_end(bucket_idx);
            let f_remainder = f.get_remainder(f_run_end);
            if let Some(b_run_end) = block.run_end(bucket_idx) {
                assert_eq!(f_run_end, b_run_end);

                if let Some(b_remainder) = block.get_remainder(b_run_end) {
                    actual_block_remainders.push(b_remainder);
                    assert_eq!(f_remainder, b_remainder);
                }
            }
        }

        assert_eq!(actual_block_remainders, expected_block_remainders);
    }

    #[test]
    fn single_block_no_next_contains() {
        let mut buckets_layout: Vec<u64> = vec![0; 10];
        buckets_layout[0] = 57;
        buckets_layout[1] = 58;
        buckets_layout[2] = 58;
        buckets_layout[3] = 59;
        buckets_layout[4] = 60;
        buckets_layout[5] = 61;
        buckets_layout[6] = 62;
        buckets_layout[7] = 63;
        buckets_layout[8] = 64;
        buckets_layout[9] = 65;

        let run_ends = [57, 59, 59, 60, 61, 62, 63];

        let (f, remainders, fingerprints) = generate_filter(100, &buckets_layout);

        let mut collision_remainders = [remainders[1], remainders[2]];
        collision_remainders.sort();
        let next = None;

        let block = extract_block(&f, 0, next);

        let mut actual_contains_count = 0;

        for fingerprint in fingerprints {
            let f_result = f.contains_fingerprint(fingerprint);
            if let Some(b_result) = block.contains_fingerprint(fingerprint) {
                actual_contains_count += 1;
                assert_eq!(f_result, b_result);
            }
        }

        assert_eq!(actual_contains_count, run_ends.len());
    }

    #[test]
    fn multiple_blocks_with_next_contains() {
        let buckets_layout = vec![55, 55, 56, 57, 57, 58, 59, 60, 60, 60, 61, 62, 62, 64, 65];

        let (f, _remainders, fingerprints) = generate_filter(100, &buckets_layout);

        let block_next = extract_block(&f, 1, None);
        let next = Some(Box::new(block_next));

        let block = extract_block(&f, 0, next);

        let mut actual_contains_count = 0;

        for &fingerprint in &fingerprints {
            let f_result = f.contains_fingerprint(fingerprint);
            if let Some(b_result) = block.contains_fingerprint(fingerprint) {
                actual_contains_count += 1;
                assert_eq!(f_result, b_result);
            }
        }

        assert_eq!(
            actual_contains_count,
            buckets_layout
                .iter()
                .filter(|&bucket_idx| *bucket_idx < 64)
                .count()
        );
    }

    #[test]
    fn extract_by_fingerprint() {
        let buckets_layout = vec![
            55, 55, 56, 57, 57, 58, 59, 60, 60, 60, 61, 62, 62, 64, 65, 127, 127,
        ];

        let (f, _remainders, fingerprints) = generate_filter(100, &buckets_layout);

        let mut fingerprints_count = 0;
        for &fingerprint in fingerprints.iter() {
            let block = Block::extract_block_by_fingerprint(&f, fingerprint);

            let f_result = f.contains_fingerprint(fingerprint);
            if let Some(b_result) = block.contains_fingerprint(fingerprint) {
                fingerprints_count += 1;
                assert_eq!(f_result, b_result);
            }
        }

        assert_eq!(fingerprints.len(), fingerprints_count);
    }
}
