#!/usr/bin/env python3
"""
Convert DDNA SNP array txt files to PLINK tped/tfam format.

Input:  Directory containing numbered sample subdirectories,
        each with a *_GSAv3-DTC_GRCh38*.txt genotype file.
Output: PLINK tped/tfam files ready for plink2 conversion.

Usage:
    python convert_ddna_to_plink.py <data_dir> <output_prefix> [min_gencall_score]
"""

import os, sys, glob

def read_ddna_file(filepath):
    """Read a DDNA genotype file, return lists of (rsid, chr, pos, genotype, gs)."""
    rsids, chroms, positions, genotypes, scores = [], [], [], [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            # Skip header line
            if parts[0] == 'rsid':
                continue
            rsids.append(parts[0])
            chroms.append(parts[1])
            positions.append(int(parts[2]))
            genotypes.append(parts[3])
            try:
                scores.append(float(parts[4]))
            except ValueError:
                scores.append(0.0)
    return rsids, chroms, positions, genotypes, scores


def main():
    data_dir = sys.argv[1]
    out_prefix = sys.argv[2]
    min_gs = float(sys.argv[3]) if len(sys.argv) > 3 else 0.15

    # Find sample directories (numeric names)
    sample_dirs = sorted(
        d for d in glob.glob(os.path.join(data_dir, '[0-9]*'))
        if os.path.isdir(d)
    )
    if not sample_dirs:
        print("ERROR: No sample directories found in", data_dir)
        sys.exit(1)

    print(f"Found {len(sample_dirs)} samples")

    # --- Build reference SNP list from first sample ---
    first_file = glob.glob(os.path.join(sample_dirs[0], '*.txt'))[0]
    print(f"Building SNP reference from: {os.path.basename(first_file)}")

    rsids, chroms, positions, genotypes, scores = read_ddna_file(first_file)

    # Filter: autosomes only, valid 2-char SNP genotypes
    valid_bases = set('ACGT')
    autosomes = set(str(c) for c in range(1, 23))
    seen_pos = set()

    keep_indices = []
    for i in range(len(rsids)):
        if chroms[i] not in autosomes:
            continue
        gt = genotypes[i]
        if len(gt) != 2 or gt[0] not in valid_bases or gt[1] not in valid_bases:
            continue
        pos_key = (chroms[i], positions[i])
        if pos_key in seen_pos:
            continue
        seen_pos.add(pos_key)
        keep_indices.append(i)

    # Sort by chromosome (numeric) then position
    keep_indices.sort(key=lambda i: (int(chroms[i]), positions[i]))

    ref_rsids = [rsids[i] for i in keep_indices]
    ref_chroms = [chroms[i] for i in keep_indices]
    ref_positions = [positions[i] for i in keep_indices]
    n_snps = len(ref_rsids)
    rsid_to_idx = {rsid: idx for idx, rsid in enumerate(ref_rsids)}

    print(f"Reference SNP set: {n_snps} autosomal SNPs")

    # --- Read all samples ---
    samples = []
    # Store alleles as byte arrays for efficiency
    alleles1 = []  # list of lists, one per sample
    alleles2 = []

    for sample_dir in sample_dirs:
        sid = os.path.basename(sample_dir)
        txt_files = glob.glob(os.path.join(sample_dir, '*.txt'))
        if not txt_files:
            print(f"  WARNING: No txt file in {sample_dir}, skipping")
            continue

        print(f"  Reading {sid}...", end=' ', flush=True)

        s_rsids, s_chroms, s_pos, s_gts, s_scores = read_ddna_file(txt_files[0])

        a1 = ['0'] * n_snps
        a2 = ['0'] * n_snps
        matched = 0

        for i in range(len(s_rsids)):
            rsid = s_rsids[i]
            if rsid not in rsid_to_idx:
                continue
            idx = rsid_to_idx[rsid]
            gt = s_gts[i]
            gs = s_scores[i]
            if gs < min_gs:
                continue
            if len(gt) != 2 or gt[0] not in valid_bases or gt[1] not in valid_bases:
                continue
            a1[idx] = gt[0]
            a2[idx] = gt[1]
            matched += 1

        pct = 100.0 * matched / n_snps
        print(f"{matched}/{n_snps} ({pct:.1f}%) genotyped")

        samples.append(sid)
        alleles1.append(a1)
        alleles2.append(a2)

    n_samples = len(samples)
    print(f"\nLoaded {n_samples} samples")

    # --- Write .tfam ---
    tfam_path = out_prefix + '.tfam'
    with open(tfam_path, 'w') as f:
        for sid in samples:
            f.write(f"0\t{sid}\t0\t0\t0\t-9\n")

    # --- Write .tped ---
    tped_path = out_prefix + '.tped'
    print(f"Writing {n_snps} SNPs x {n_samples} samples to tped...")
    with open(tped_path, 'w') as f:
        for i in range(n_snps):
            parts = [ref_chroms[i], ref_rsids[i], '0', str(ref_positions[i])]
            for j in range(n_samples):
                parts.append(alleles1[j][i])
                parts.append(alleles2[j][i])
            f.write('\t'.join(parts) + '\n')

            if (i + 1) % 100000 == 0:
                print(f"  {i+1}/{n_snps} SNPs written...")

    print(f"\nDone!")
    print(f"  {tfam_path} ({n_samples} samples)")
    print(f"  {tped_path} ({n_snps} SNPs)")


if __name__ == '__main__':
    main()
