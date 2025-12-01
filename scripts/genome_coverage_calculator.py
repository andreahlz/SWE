genomes = {
    "A_muciniphila": {"size": 2900000, "abundance": 0.08},
    "B_fragilis": {"size": 5200000, "abundance": 0.40},
    "B_longum": {"size": 2400000 , "abundance": 0.05},
    "B_obeum": {"size": 3600000, "abundance": 0.15},
    "E_faecalis": {"size": 2900000, "abundance": 0.06},
    "F_prausnitzii": {"size": 3100000, "abundance": 0.20},
    "P_intermedia": {"size": 2700000, "abundance": 0.04},
    "R_intestinalis": {"size": 4500000, "abundance": 0.02}
}

target_coverage = 10
short_read_length = 150
long_read_length = 10000

print("=" * 60)
print("METAGENOMIC READ CALCULATION FOR 10X COVERAGE")
print("=" * 60)

print("\n### SHORT READS (150 bp) ###")
total_short_reads = 0
for species, data in genomes.items():
    effective_size = data["size"]*data["abundance"]
    reads_needed = (target_coverage * effective_size) / short_read_length
    total_short_reads += reads_needed
    print(f"{species:40} {int(reads_needed):>10,} reads")


print("\n### LONG READS (10 kb) ###")
total_long_reads = 0
for species, data in genomes.items():
    effective_size = data["size"]*data["abundance"]
    reads_needed = (target_coverage * effective_size) / long_read_length
    total_long_reads += reads_needed
    print(f"{species:40} {int(reads_needed):>10,} reads")

print(f"\n{'TOTAL LONG READS':40} {int(total_long_reads):>10,} reads")