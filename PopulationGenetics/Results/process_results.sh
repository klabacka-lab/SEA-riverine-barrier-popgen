#!/usr/bin/env bash
set -euo pipefail

out_csv="fst_summary.csv"

echo "Genus,specific_epithet,gene,river,FST,greater_perm,total_perm,p-value" > "$out_csv"

for file in Output_*.txt; do
    # Skip if no matching files
    [[ -e "$file" ]] || continue

    # Remove prefix and suffix, then split on underscores
    base="${file#Output_}"
    base="${base%.txt}"

    IFS='_' read -r genus specific_epithet gene river <<< "$base"

    # Extract values from the first 4 lines
    fst=$(awk -F': ' '/^True FST:/ {print $2; exit}' "$file")
    greater_perm=$(awk -F': ' '/^# permutations with FST > true FST:/ {print $2; exit}' "$file")
    total_perm=$(awk -F': ' '/^Total permutations:/ {print $2; exit}' "$file")
    p_value=$(awk -F': ' '/^p-value:/ {print $2; exit}' "$file")

    # Write CSV row
    printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$genus" "$specific_epithet" "$gene" "$river" \
        "$fst" "$greater_perm" "$total_perm" "$p_value" >> "$out_csv"
done

echo "Wrote $out_csv"
