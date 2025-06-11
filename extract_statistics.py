import re
import sys
import subprocess
import io
import csv
import os
import glob

# get number of raw reads from fastqc
def fastqc_read_count(filename):
    if not filename or not os.path.exists(filename):
        return "NA"

    with open(filename, 'r') as f:

        for line in f:

            match = re.search(r'Total Sequences</td><td>([\d,]+)', line)

            if match:

                return int(match.group(1).replace(',', ''))

    raise ValueError(f"Raw read count not found in {filename}")


# get number of mapped reads from flagstat file
def get_mapped_reads(filename):
    if not filename or not os.path.exists(filename):
        return "NA"

    with open(filename, 'r') as f:

        for line in f:

            match = re.match(r'^(\d+)\s+\+\s+\d+\s+in total', line)

            if match:

                return int(match.group(1))

    raise ValueError(f"Mapped read count not found in {filename}")


# get the read length from AdapterRemoval settings file
def get_read_length(filename):
    if not filename or not os.path.exists(filename):
        return "NA"

    with open(filename, 'r') as f:

        for line in f:

            match = re.match(r'^Average read length of trimmed reads:\s+([0-9.]+)', line)

            if match:

                return float(match.group(1))

    raise ValueError(f"Read length not found in {filename}")



# samtools coverage & statistics

def get_coverage_stats(bam_file, mapped_reads):
    if not bam_file or not os.path.exists(bam_file):
        return {
            'genome_coverage_depth': "NA",
            'genome_coverage_breadth': "NA",
            'mtdna_reads': "NA",
            'depth_mtdna': "NA",
            'breadth_mtdna': "NA",
            'mt_nuc_ratio': "NA"
        }

    if mapped_reads == "NA":
        return {
            'genome_coverage_depth': "NA",
            'genome_coverage_breadth': "NA",
            'mtdna_reads': "NA",
            'depth_mtdna': "NA",
            'breadth_mtdna': "NA",
            'mt_nuc_ratio': "NA"
        }

    try:
        result = subprocess.run(
            ["samtools-1.18-rocky", "coverage", bam_file],
            check=True,
            text=True,
            stdout=subprocess.PIPE
        )

        total_pos = 0
        total_bases_cov = 0
        total_bases_aligned = 0

        mtdna_reads = 0
        mtdna_depth = 0.0
        mtdna_breadth = 0.0

        for line in result.stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split('\t')
            rname = fields[0]

            endposition = int(fields[2])
            numreads = int(fields[3])
            covbases = int(fields[4])
            coverage = float(fields[5])
            meandepth = float(fields[6])

            total_pos += endposition
            total_bases_cov += covbases
            total_bases_aligned += endposition * meandepth

            if rname == "chrM":
                mtdna_reads = numreads
                mtdna_depth = meandepth
                mtdna_breadth = coverage

        if total_pos == 0:
            raise ValueError("No data rows found or total endpos is zero.")

        coverage_depth = total_bases_aligned / total_pos
        coverage_breadth = (total_bases_cov / total_pos) * 100

        if mapped_reads - mtdna_reads > 0:
            mt_nuc_ratio = mtdna_reads / (mapped_reads - mtdna_reads)
        else:
            mt_nuc_ratio = 0.0

        return {
            'genome_coverage_depth': coverage_depth,
            'genome_coverage_breadth': coverage_breadth,
            'mtdna_reads': mtdna_reads,
            'depth_mtdna': mtdna_depth,
            'breadth_mtdna': mtdna_breadth,
            'mt_nuc_ratio': mt_nuc_ratio
        }

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"samtools coverage failed: {e}")
    except Exception as e:
        raise RuntimeError(f"Failed to parse samtools coverage output: {e}")




#get the insert size from "samtools stats"
def get_insert_size(bam_file):
    if not bam_file or not os.path.exists(bam_file):
        return "NA"

    try:
        # Compose the shell command
        cmd = (
            f"samtools-1.18-rocky stats {bam_file} | "
            "grep ^IS | cut -f 3- | "
            "awk '{sum += $1 * $2; total += $2} END {if (total > 0) print sum / total; else print \"0\"}'"
        )

        # Run it through shell and capture the output
        result = subprocess.run(cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE)

        return float(result.stdout.strip())

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Insert size command failed: {e}")
    except ValueError:
        raise ValueError(f"Could not parse insert size from output: {result.stdout}")


#get number of duplicated reads from flagstat file
def get_duplicated_reads(filename):
    if not filename or not os.path.exists(filename):
        return "NA"

    with open(filename, 'r') as f:

        for line in f:

            match = re.match(r'^(\d+)\s+\+\s+\d+\s+duplicates', line)

            if match:

                return int(match.group(1))

    raise ValueError(f"Duplicated read count not found in {filename}")







#######APPLY THE FUNCTIONS######

def process_samples(base_dir, sample_list_file):

#ADD COLUMNS HERE
    print("\t".join([

        "Sample", "Reads Lane 1", "Reads Lane 2", "Total Raw Reads", "Mapped Reads", "% Mapped", "Discarded Reads", "Average Length",

        "Depth of Coverage", "Breadth of Coverage", "mtDNA reads", "Depth of Coverage for mtDNA", "Breadth of Coverage for mtDNA", 

        "Mitochondrial to Nuclear Read Ratio", "Mean Insert Size", "Duplicated Reads", "Duplication Rate"
    ]))


    with open(sample_list_file, 'r') as f:

        sample_ids = [line.strip() for line in f if line.strip()]

        for sample_id in sample_ids:

            sample_dir = os.path.join(base_dir, sample_id)


            def first_match(pattern, description):
                matches = glob.glob(pattern)
                if not matches:
                    return None  # Let downstream handle it as "NA"
                return matches[0]

            fastqc_L001_html = first_match(os.path.join(sample_dir, f"{sample_id}_S*_L001_R1_001_fastqc.html"), "FastQC L001 HTML")
            fastqc_L002_html = first_match(os.path.join(sample_dir, f"{sample_id}_S*_L002_R1_001_fastqc.html"), "FastQC L002 HTML")
            adapter_removal_L001_settings = first_match(os.path.join(sample_dir, f"{sample_id}_S*_L001_R1.settings"), "AdapterRemoval L001 settings")
            adapter_removal_L002_settings = first_match(os.path.join(sample_dir, f"{sample_id}_S*_L002_R1.settings"), "AdapterRemoval L002 settings")

            markdup_flagstat_file = os.path.join(sample_dir, f"{sample_id}.merged.sorted.markDup.statistics")
            bam_path = os.path.join(sample_dir, f"{sample_id}.merged.sorted.markDup.sorted.recalibrated.bam")


            stats = {
                'lane1_reads': fastqc_read_count(fastqc_L001_html),
                'lane2_reads': fastqc_read_count(fastqc_L002_html)
            }

            # Ensure total_reads propagates NA correctly
            if stats['lane1_reads'] == "NA" or stats['lane2_reads'] == "NA":
                stats['total_reads'] = "NA"
            else:
                stats['total_reads'] = 2 * stats['lane1_reads'] + 2 * stats['lane2_reads']


            stats['mapped_reads'] = get_mapped_reads(markdup_flagstat_file)

            # Handle mapped reads and dependent calculations
            if stats['total_reads'] == "NA":
                stats['percent_mapped'] = "NA"
                stats['discarded_reads'] = "NA"
            else:
                stats['percent_mapped'] = (stats['mapped_reads'] / stats['total_reads'] * 100) if stats['total_reads'] != 0 else "NA"
                stats['discarded_reads'] = stats['total_reads'] - stats['mapped_reads']

            stats['average_length'] = (get_read_length(adapter_removal_L001_settings) + get_read_length(adapter_removal_L002_settings))/2

            coverage_stats = get_coverage_stats(bam_path, stats['mapped_reads'])
            stats.update(coverage_stats)

            stats['mean_insert_size'] = get_insert_size(bam_path)
            stats['duplicated_reads'] = get_duplicated_reads(markdup_flagstat_file)
            stats['duplication_rate'] = stats['duplicated_reads']/stats['mapped_reads'] * 100


            print("\t".join(map(str, [

                sample_id,

                stats['lane1_reads'],
                
                stats['lane2_reads'],
                
                stats['total_reads'],
                
                stats['mapped_reads'],
                
                f"{stats['percent_mapped']:.2f}" if isinstance(stats['percent_mapped'], (int, float)) else "NA",
                
                stats['discarded_reads'],

                f"{stats['average_length']:.2f}" if isinstance(stats['average_length'], (int, float)) else "NA",
                
                f"{stats['genome_coverage_depth']:.2f}" if isinstance(stats['genome_coverage_depth'], (int, float)) else "NA",
                
                f"{stats['genome_coverage_breadth']:.2f}" if isinstance(stats['genome_coverage_breadth'], (int, float)) else "NA",
                
                stats['mtdna_reads'],
                
                f"{stats['depth_mtdna']:.2f}" if isinstance(stats['depth_mtdna'], (int, float)) else "NA",
                
                f"{stats['breadth_mtdna']:.2f}" if isinstance(stats['breadth_mtdna'], (int, float)) else "NA",
                
                f"{stats['mt_nuc_ratio']:.4f}" if isinstance(stats['mt_nuc_ratio'], (int, float)) else "NA",
                
                f"{stats['mean_insert_size']:.2f}" if isinstance(stats['mean_insert_size'], (int, float)) else "NA",
                
                stats['duplicated_reads'],
                
                f"{stats['duplication_rate']:.2f}" if isinstance(stats['duplication_rate'], (int, float)) else "NA",


        ])))




if __name__ == "__main__":

    if len(sys.argv) != 3:

        print("Usage: python summarize_qualimap.py <base_dir> <sample_list.txt>")

        sys.exit(1)

    base_dir = sys.argv[1]
    sample_list_file = sys.argv[2]

    process_samples(base_dir, sample_list_file)
