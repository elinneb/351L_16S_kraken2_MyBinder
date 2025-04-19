BEGIN {FS=OFS="\t"}
NR == 1 {
    print $0  # Print header unchanged
    next
}
{
    # Calculate row sum (skip first column)
    sum = 0
    for (i=2; i<=NF; i++) sum += $i
    
    # Store row and sum
    rows[NR] = $0
    sums[NR] = sum
    original_order[NR] = NR
}
END {
    # Simple bubble sort by sum (descending)
    for (i=2; i<=NR; i++) {
        for (j=i+1; j<=NR; j++) {
            if (sums[i] < sums[j]) {
                # Swap sums
                tmp = sums[i]
                sums[i] = sums[j]
                sums[j] = tmp
                # Swap row indices
                tmp = original_order[i]
                original_order[i] = original_order[j]
                original_order[j] = tmp
            }
        }
    }
    
    # Print sorted rows
    for (i=2; i<=NR; i++) {
        print rows[original_order[i]]
    }
}
