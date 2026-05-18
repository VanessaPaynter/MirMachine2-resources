#!/usr/bin/env bash
# Usage: ./filter_dotplots_perfectpair.sh input.dotplot output_nonperfect.dotplot
# - Prints ALL records with debug info
# - Always checks structure length and bracket continuity
# - Writes ONLY non-perfect-pairing records to output file

input="$1"
output="$2"

awk -v out="$output" '
BEGIN {
    RS=">"; ORS="";
    if (out == "") {
        print "Usage: ./filter_dotplot_perfectpairing_struct.sh input.dotplot output.dotplot\n" > "/dev/stderr";
        exit 1;
    }
}
NR > 1 {
    n = split($0, lines, "\n")

    header  = ">" lines[1]
    seq     = lines[2]
    struct  = lines[3]

    gsub(/[[:space:]]+$/, "", struct)
	# Extract only structure before first space (to remove free energy part)
	split(struct, parts, " ")
	struct_clean = parts[1]
	struct_len = length(struct_clean)

    print "---- Checking record ----\n" > "/dev/stderr"
    print header > "/dev/stderr"
    print "   Structure length: " struct_len ((struct_len > 100) ? "  [>100]" : "") > "/dev/stderr"

    # Check for long continuous bracket region
	# First: original fast regex check
	regex_hit = (struct_clean ~ /\({15,}/ || struct_clean ~ /\){15,}/)

	has_perfect = 0

	if (regex_hit) {

	    # --- Parse base pairs with stack ---
	    delete pair
	    delete stack
	    top = 0

	    for (i = 1; i <= struct_len; i++) {
	        c = substr(struct_clean, i, 1)

	        if (c == "(") {
	            stack[++top] = i
	        }
	        else if (c == ")") {
	            if (top > 0) {
	                j = stack[top--]
	                pair[i] = j
	                pair[j] = i
	            }
	        }
	    }

	    # --- Check for true symmetric 15bp run ---
	    for (i = 1; i <= struct_len - 14; i++) {
	        if (i in pair) {
	            run = 1
	            left = i
	            right = pair[i]

	            for (k = 1; k < 15; k++) {
	                if ((left + k) in pair &&
	                    pair[left + k] == right - k) {
	                    run++
	                } else {
	                    break
	                }
	            }

	            if (run >= 15) {
	                has_perfect = 1
	                break
	            }
	        }
	    }
	}				

	if (regex_hit && has_perfect) {

	    print "→ PERFECT pairing detected (≥15 consecutive brackets)" > "/dev/stderr"
	    print "→ Structural integrity confirmed (true ≥15bp symmetric stem)\n" > "/dev/stderr"
	    print header "\n" seq "\n" struct "\n" > "/dev/stderr"

	}
	else if (regex_hit && !has_perfect) {

	    print "→ ≥15 consecutive brackets detected (regex hit)" > "/dev/stderr"
	    print "→ Structural integrity FAILED (not symmetric 15bp stem; likely nested/bulged)\n" > "/dev/stderr"
	    print header "\n" seq "\n" struct "\n" >> out

	}
	else {

	    print "→ No perfect pairing (kept)\n" > "/dev/stderr"
	    print header "\n" seq "\n" struct "\n" >> out
	}

	print "" > "/dev/stderr"
	}

' "$input"
