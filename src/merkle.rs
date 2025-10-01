use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// Compute SHA-256 of some data and return it as a `[u8; 32]`.
#[inline(always)]
pub fn sha256(data: &[u8]) -> [u8; 32] {
    Sha256::digest(data).into()
}

/// Build per-byte leaf hashes for a file.  Each byte is hashed on its own so
/// that inclusion proofs can be made for arbitrary slice lengths.
#[inline]
pub fn build_byte_leaf_hashes(file_bytes: &[u8]) -> Vec<[u8; 32]> {
    file_bytes.iter().map(|b| sha256(&[*b])).collect()
}

#[inline]
pub fn leaf_hashes_from_iter<I, T>(items: I) -> Vec<[u8; 32]>
where
    I: Iterator<Item = T>,
    T: AsRef<[u8]>,
{
    items.map(|item| sha256(item.as_ref())).collect()
}

#[derive(Deserialize, Serialize)]
pub struct MerkleRangeProof {
    pub left: Vec<[u8; 32]>,
    pub right: Vec<[u8; 32]>,
    pub range: (usize, usize),
}

impl MerkleRangeProof {
    pub fn verify(
        &self,
        leaf_hashes: &[[u8; 32]],
        left_index: usize,
        right_index: usize,
    ) -> [u8; 32] {
        let mut left_index_level = left_index;
        let mut right_index_level = right_index;

        let proof_path_len = self.right.len().max(self.left.len());
        let empty_witness = [0u8; 32];

        let mut left_witness = self.left.clone();
        let mut right_witness = self.right.clone();

        let mut current_level = leaf_hashes.to_vec();

        for i in 0..proof_path_len {
            let left = left_witness[i];
            let right = right_witness[i];
            let mut next_level = Vec::new();

            if left != empty_witness {
                current_level.insert(0, left);
            }
            if right != empty_witness {
                current_level.push(right);
            }

            for j in (0..current_level.len()).step_by(2) {
                let left_witness = current_level[j];
                let right_witness = current_level[j + 1];

                let combined = [&left_witness[..], &right_witness[..]].concat();
                let mut hasher = Sha256::new();
                hasher.update(&combined);
                next_level.push(hasher.finalize().into());
            }

            left_index_level /= 2;
            right_index_level /= 2;
            current_level = next_level;
        }

        current_level.last().cloned().unwrap()
    }

    pub fn build(
        leaf_hashes: &[[u8; 32]],
        start_index: usize,
        end_index: usize,
    ) -> MerkleRangeProof {
        let mut left: Vec<[u8; 32]> = Vec::new();
        let mut right: Vec<[u8; 32]> = Vec::new();

        let mut current_level = leaf_hashes.to_vec();

        let mut left_range_bound = start_index;
        let mut right_range_bound = end_index - 1;

        while current_level.len() > 1 {
            let mut next_level = Vec::new();

            if current_level.len() % 2 == 1 && !current_level.is_empty() {
                let last = current_level.last().cloned().unwrap();
                current_level.push(last);
            }

            for i in (0..current_level.len()).step_by(2) {
                let left_witness = current_level[i];
                let right_witness = current_level[i + 1];

                if i + 1 >= left_range_bound {
                    let empty_witness = [0u8; 32];
                    // [- -] [- -] [- +] [+ +] [- -]
                    //             [^ ^]
                    //             [i; i + 1]
                    if i < left_range_bound {
                        left.push(current_level[i]);
                    // [- -] [- -] [+ +] [+ +] [- -]
                    //             [^ ^]
                    //             [i; i + 1]
                    } else if i == left_range_bound {
                        // [0 1] [2 3] [4 5] [6 7] [8 9]
                        // [- -] [- -] [+ +] [+ +] [- -]
                        //             [^ ^]
                        //             [i; i + 1]
                        // if i % 2 != 1 {
                        left.push(empty_witness);
                        // [0 1] [2 3] [4 5] [6 7] [8 9]
                        // [- -] [- -] [+ +] [+ -] [- -]
                        //             [^   ^]
                        //         [i; i + 1]
                        // } else {
                        //     self.left.push(left_witness);
                        // }
                    }
                }

                if i <= right_range_bound {
                    let empty_witness = [0u8; 32];
                    // [- -] [- -] [+ +] [+ -] [- -]
                    //                   [^ ^]
                    //                   [i; i + 1]
                    if i + 1 > right_range_bound {
                        right.push(current_level[i + 1]);

                    // [- -] [- -] [+ +] [+ +] [- -]
                    //                   [^ ^]
                    //                   [i; i + 1]
                    } else if i + 1 == right_range_bound {
                        // [0 1] [2 3] [4 5] [6 7] [8 9]
                        // [- -] [- -] [+ +] [+ +] [- -]
                        //                   [^ ^]
                        //                   [i; i + 1]
                        // if i % 2 != 1 {
                        right.push(empty_witness);
                        // [0 1] [2 3] [4 5] [6 7] [8 9]
                        // [- -] [- -] [+ +] [+ -] [- -]
                        //               [^   ^]
                        //               [i; i + 1]
                        // } else {
                        //     self.right.push(right_witness);
                        // }
                    }
                }

                let combined = [&left_witness[..], &right_witness[..]].concat();
                let mut hasher = Sha256::new();
                hasher.update(&combined);
                next_level.push(hasher.finalize().into());
            }

            left_range_bound /= 2;
            right_range_bound /= 2;
            current_level = next_level;
        }

        MerkleRangeProof {
            left,
            right,
            range: (start_index, end_index),
        }
    }
}

/// Compute the Merkle root from a vector of leaf hashes.  The tree uses
/// SHA-256 and left-concatenation ordering.  If the number of nodes at any
/// level is odd, the last node is duplicated (standard Bitcoin/ETH style).
#[inline]
pub fn compute_root(mut level: Vec<[u8; 32]>) -> [u8; 32] {
    if level.is_empty() {
        return [0u8; 32];
    }
    if level.len() == 1 {
        return level[0];
    }

    while level.len() > 1 {
        if level.len() % 2 == 1 {
            // Duplicate last node to make count even.
            let last = *level.last().unwrap();
            level.push(last);
        }
        let mut next_level = Vec::with_capacity((level.len() + 1) / 2);
        for pair in level.chunks_exact(2) {
            let combined = [&pair[0][..], &pair[1][..]].concat();
            next_level.push(sha256(&combined));
        }
        level = next_level;
    }
    level[0]
}

/// Build the sibling path needed to prove inclusion of `leaf_index`.
#[inline]
pub fn build_proof_path(leaves: &[[u8; 32]], mut leaf_index: usize) -> Vec<[u8; 32]> {
    if leaves.len() <= 1 {
        return Vec::new();
    }

    let mut proof = Vec::new();
    let mut current_level: Vec<[u8; 32]> = leaves.to_vec();

    while current_level.len() > 1 {
        if current_level.len() % 2 == 1 {
            let last = *current_level.last().unwrap();
            current_level.push(last);
        }

        let mut next_level = Vec::with_capacity((current_level.len() + 1) / 2);
        for chunk in current_level.chunks_exact(2).enumerate() {
            let (pair_index, pair) = chunk;
            let left = pair[0];
            let right = pair[1];

            let global_left_index = pair_index * 2;
            if leaf_index == global_left_index {
                // Current node is the left child.
                proof.push(right);
            } else if leaf_index == global_left_index + 1 {
                // Current node is the right child.
                proof.push(left);
            }

            let combined = [&left[..], &right[..]].concat();
            next_level.push(sha256(&combined));
        }
        leaf_index /= 2;
        current_level = next_level;
    }

    proof
}

/// Convenience wrapper that returns `(root, proof_path)` for a single leaf.
#[inline]
pub fn root_and_path_for_leaf(leaves: &[[u8; 32]], leaf_index: usize) -> ([u8; 32], Vec<[u8; 32]>) {
    let root = compute_root(leaves.to_vec());
    let path = build_proof_path(leaves, leaf_index);
    (root, path)
}

/// Build a complete Merkle tree from a byte slice using 32-byte chunks as
/// leaves and return the Merkle root together with the vector of leaf hashes.
///
/// The algorithm guarantees that the same input bytes always produce the same
/// root hash.  Chunks shorter than 32 bytes (only possible for the final
/// chunk) are zero-padded on the right.
///
/// This helper exists in addition to [`build_byte_leaf_hashes`] which operates
/// on *per-byte* leaves and is useful when proofs need to be made for arbitrary
/// slice lengths.  `build_merkle_tree` is more efficient for use-cases that can
/// work with 32-byte alignment (e.g. file-level proofs).
#[inline]
pub fn build_merkle_tree(file_bytes: &[u8]) -> ([u8; 32], Vec<[u8; 32]>) {
    // 1. Build the leaf hashes ------------------------------------------------
    let leaves: Vec<[u8; 32]> = file_bytes
        .chunks(32)
        .map(|chunk| {
            let mut padded = [0u8; 32];
            padded[..chunk.len()].copy_from_slice(chunk);
            sha256(&padded)
        })
        .collect();

    // Single-leaf special-case -------------------------------------------------
    if leaves.len() == 1 {
        return (leaves[0], leaves);
    }

    // 2. Build the Merkle tree bottom-up --------------------------------------
    let mut current_level = leaves.clone();
    while current_level.len() > 1 {
        // Duplicate the last node to ensure an even number of siblings.
        if current_level.len() % 2 == 1 {
            current_level.push(*current_level.last().unwrap());
        }

        let mut next_level = Vec::with_capacity((current_level.len() + 1) / 2);
        for pair in current_level.chunks_exact(2) {
            let combined = [&pair[0][..], &pair[1][..]].concat();
            next_level.push(sha256(&combined));
        }
        current_level = next_level;
    }

    (current_level[0], leaves)
}
