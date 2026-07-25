#!/bin/bash

template="s.run_popgen_analysis.sh"
list="alignment_list.txt"

while IFS= read -r name; do
    # Skip empty lines
    [[ -z "$name" ]] && continue

    outfile="s.run_popgen_analysis_${name}.sh"

    # Replace the placeholder and write the new script
    sed "s/<replace_me>/${name}/g" "$template" > "$outfile"

    # Make the new script executable (optional)
    chmod +x "$outfile"

    echo "Created $outfile"
done < "$list"
