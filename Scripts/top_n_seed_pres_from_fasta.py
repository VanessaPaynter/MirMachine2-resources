import pandas as pd
import re
import os
from glob import glob

###################################
# USER INPUTS
###################################

DEBUG = True

input_directory = "/Users/vpa021/Downloads/MIR-12_rerun_May26/MIR-12_fastas"

# CSV outputs
output_all_csv = "seed_plus_nonseed_Mir12_lep_Jackdata_MM2_topN.csv"
output_seeds_only_csv = "seeds_only_Mir12_lep_Jackdata_MM2_topN.csv"
output_no_seed_csv = "nonseeds_only_no_seed_families_Mir12_lep_Jackdata_MM2.csv"

# FASTA outputs
output_all_fasta = "seed_plus_nonseed_Mir12_lep_Jackdata_MM2_topN.fasta"
output_seeds_only_fasta = "seeds_only_Mir12_lep_Jackdata_MM2_topN.fasta"
output_no_seed_fasta = "nonseeds_only_no_seed_families_Mir12_lep_Jackdata_MM2.fasta"

# paralogue numbers per family
paralogue_numbers = {
    'Bantam': 1,
    'Iab-4': 2,
    'Let-7': 1,
    'Mir-1': 1,
    'Mir-10': 4,
    'Mir-1000': 1,
    'Mir-1001': 1,
    'Mir-1003': 1,
    'Mir-1006': 1,
    'Mir-1007': 1,
    'Mir-1010': 1,
    'Mir-10485': 1,
    'Mir-1175': 1.5,
    'Mir-12': 1,
    'Mir-124': 1,
    'Mir-133': 1,
    'Mir-137': 1,
    'Mir-14': 1,
    'Mir-184': 1,
    'Mir-190': 1.5,
    'Mir-193': 3.5,
    'Mir-2': 9,
    'Mir-2001': 1,
    'Mir-210': 2,
    'Mir-216': 2,
    'Mir-219': 2,
    'Mir-22': 1,
    'Mir-2279': 1,
    'Mir-2499': 1,
    'Mir-2501': 1,
    'Mir-252': 3.5,
    'Mir-275': 1,
    'Mir-2755': 1,
    'Mir-2756': 1,
    'Mir-276': 2,
    'Mir-2763': 2,
    'Mir-2765': 1,
    'Mir-2766': 2,
    'Mir-2767': 1,
    'Mir-277': 1,
    'Mir-278': 1,
    'Mir-2786': 1,
    'Mir-279': 3.5,
    'Mir-2796': 1,
    'Mir-281': 3.5,
    'Mir-282': 1,
    'Mir-284': 2,
    'Mir-29': 2.5,
    'Mir-3': 9,
    'Mir-303': 1,
    'Mir-305': 1,
    'Mir-306': 1,
    'Mir-31': 1.5,
    'Mir-314': 1,
    'Mir-315': 1,
    'Mir-316': 1,
    'Mir-317': 1,
    'Mir-33': 1.5,
    'Mir-3327': 1,
    'Mir-3338': 1,
    'Mir-34': 1.5,
    'Mir-375': 1,
    'Mir-4969': 1,
    'Mir-4983': 1,
    'Mir-6301': 2,
    'Mir-6302': 2,
    'Mir-6304': 1,
    'Mir-6305': 1,
    'Mir-6307': 1,
    'Mir-67': 2.5,
    'Mir-7': 1,
    'Mir-71': 1,
    'Mir-750': 1,
    'Mir-76': 1,
    'Mir-8': 1,
    'Mir-87': 1.5,
    'Mir-9': 2.5,
    'Mir-90': 1,
    'Mir-91': 1,
    'Mir-92': 5,
    'Mir-927': 2,
    'Mir-929': 1,
    'Mir-93': 1,
    'Mir-932': 1,
    'Mir-9388': 1,
    'Mir-955': 1,
    'Mir-956': 2,
    'Mir-957': 1,
    'Mir-959': 2,
    'Mir-96': 2,
    'Mir-960': 1,
    'Mir-961': 1,
    'Mir-962': 1,
    'Mir-963': 2,
    'Mir-964': 2,
    'Mir-965': 1.5,
    'Mir-966': 1,
    'Mir-967': 1,
    'Mir-969': 2,
    'Mir-970': 1,
    'Mir-971': 1.5,
    'Mir-972': 1,
    'Mir-973': 1,
    'Mir-974': 2,
    'Mir-975': 1,
    'Mir-976': 1,
    'Mir-977': 1,
    'Mir-978': 2,
    'Mir-982': 1,
    'Mir-983': 2,
    'Mir-984': 1,
    'Mir-985': 1,
    'Mir-986': 1,
    'Mir-987': 1,
    'Mir-988': 1,
    'Mir-989': 1.5,
    'Mir-991': 1,
    'Mir-992': 2,
    'Mir-997': 1,
    'Mir-999': 1
}

###################################
# HEADER PARSER
###################################

def parse_header(header):
    family, rest = header.split(".", 1)

    pattern = (
        r"PRE_(.+?)_"          # chromosome/contig (non-greedy)
        r"(\d+)_"           # start
        r"(\d+)_"           # end
        r"\(([-+])\)_"      # strand
        r"([\d.]+)_"        # bitscore
        r"(HIGHconf|LOWconf)"  # confidence
    )

    match = re.search(pattern, rest)
    if not match:
        print("FAILED HEADER:", header)
        raise ValueError(f"Header parsing failed: {header}")

    chromosome = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    strand = match.group(4)
    bitscore = float(match.group(5))
    confidence = match.group(6)

    # Extract seeds
    seed_pattern = r"seed\((.*?)\)"
    seed_matches = re.findall(seed_pattern, header)
    seeds = []
    for m in seed_matches:
        seeds.extend(m.split(","))
    seeds = [s.strip() for s in seeds if s.strip()]
    seed_present = len(seeds) > 0

    return {
        "Family": family,
        "Chromosome": chromosome,
        "Start": start,
        "End": end,
        "Strand": strand,
        "Bitscore": bitscore,
        "Confidence": confidence,
        "SeedPresent": seed_present,
        "Seeds": ",".join(seeds) if seeds else "",
        "Header": header
    }

###################################
# LOAD FASTA FILES
###################################

records = []
fasta_files = glob(os.path.join(input_directory, "*.fasta"))

for fasta_file in fasta_files:
    gca = os.path.splitext(os.path.basename(fasta_file))[0]
    with open(fasta_file, "r") as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                sequence = line
                parsed = parse_header(header)
                parsed["Sequence"] = sequence
                parsed["GCA"] = gca
                records.append(parsed)

data = pd.DataFrame(records)

###################################
# SYMMETRIC TOP-N LOGIC WITH DEBUG MSGs
###################################

all_topN = []
seeds_only = []
no_seed_nonseeds = []

for (family, gca), group in data.groupby(["Family", "GCA"]):

    if family not in paralogue_numbers:
        continue

    N = round(paralogue_numbers[family])
    seeds = group[group["SeedPresent"]].sort_values("Bitscore", ascending=False).reset_index(drop=True)
    nonseeds = group[~group["SeedPresent"]].sort_values("Bitscore", ascending=False).reset_index(drop=True)

    if DEBUG:
        print("\n=======================================")
        print(f"Family: {family} | Genome: {gca} | Expected N={N}")

    if not seeds.empty:
        if DEBUG:
            print("Seeds detected (sorted by bitscore):")
            print(seeds[["Bitscore", "Seeds"]])

        retained = seeds.head(N).copy()

        if len(retained) < N:
            weakest_seed_score = retained["Bitscore"].min()
            needed = N - len(retained)
            eligible_nonseeds = nonseeds[nonseeds["Bitscore"] > weakest_seed_score].reset_index(drop=True)

            if DEBUG:
                print(f"Weakest retained seed score: {weakest_seed_score}")
                print(f"Need {needed} additional entries from non-seeds")
                print("All non-seeds (sorted by bitscore):")
                print(nonseeds[["Bitscore", "Seeds"]])
                print("Eligible non-seeds (must be > weakest seed):")
                print(eligible_nonseeds[["Bitscore", "Seeds"]])

            retained = pd.concat([retained, eligible_nonseeds.head(needed)]).reset_index(drop=True)

        if DEBUG:
            print("Final retained top-N entries for this group:")
            print(retained[["Bitscore", "SeedPresent", "Seeds"]])

        all_topN.append(retained)
        seeds_only.append(retained[retained["SeedPresent"]].reset_index(drop=True))

    else:
        if DEBUG:
            print("No seeds detected.")
            print("All non-seeds (sorted by bitscore):")
            print(nonseeds[["Bitscore", "Seeds"]])

        retained = nonseeds.head(N).copy().reset_index(drop=True)
        if DEBUG:
            print(f"Top {len(retained)} non-seeds retained:")
            print(retained[["Bitscore", "SeedPresent", "Seeds"]])

        all_topN.append(retained)
        no_seed_nonseeds.append(retained)

###################################
# SAVE CSVs
###################################

df_all = pd.concat(all_topN, ignore_index=True)
df_seeds_only = pd.concat(seeds_only, ignore_index=True) if seeds_only else pd.DataFrame()
df_no_seed = pd.concat(no_seed_nonseeds, ignore_index=True) if no_seed_nonseeds else pd.DataFrame()

df_all.to_csv(output_all_csv, index=False)
df_seeds_only.to_csv(output_seeds_only_csv, index=False)
df_no_seed.to_csv(output_no_seed_csv, index=False)

###################################
# SAVE FASTA FILES
###################################

def save_fasta(df, filename):
    if df.empty:
        return
    with open(filename, "w") as f:
        for _, row in df.iterrows():
            # Prepend GCA/genome name to header
            f.write(f">{row['GCA']}|{row['Header']}\n{row['Sequence']}\n")

save_fasta(df_all, output_all_fasta)
save_fasta(df_seeds_only, output_seeds_only_fasta)
save_fasta(df_no_seed, output_no_seed_fasta)

print("\nFinished processing.")
print(f"CSV & FASTA files written for all three categories.")
print("Total topN retained:", len(df_all))
print("Seeds only:", len(df_seeds_only))
print("Non-seeds (no-seed families):", len(df_no_seed))