"""
rate_partitions.py
Written by: T Malm - 20120922
Script to read a TIGER rate file (or other rate files) and produce partitions
from the rates usable for analyses.

Currently partition sizes are calculated as follows:
The first partition is calculated as

highest_rate - ((highest_rate - minimum_rate)/(divfactor)), the remaining bins as:
higher_rate - ((higher_rate - lower_rate)/divfactor + partition_number * 0.3)

leading to smaller bin rate spans with faster rates (=lower rates), the last
partition will be created when 10% or less characters are left.

Run the script from terminal with a rate-file, as first argument, and the wanted
divfactor, as second argument, and it will produce a new output file with
partition rate summaries and partition schemes for MrBayes and PHYLIP.

The rate file should be a file with rates separated by hard return (TIGER rate
output -rl command) and in order of the character alignment you want to use later.

Run as: "python rate_partitions.py X Y" where X is the input rate file and Y
is the division factor (e.g. 1.35)

The output file will be named after the input values: inputfilename_divfactor.txt
"""
import argparse
import os
import re
import sys


def run(infile, divnum):
    """Executes computations and returns contents to be written in output file"""
    input_data = read_input_file(infile)
    infile_basename = os.path.basename(infile)
    
    max_rate = max(input_data)
    min_rate = min(input_data)
    num_chars = len(input_data)
    print("\nTotal data - {} sites.\n".format(num_chars))
    print("Slowest rate: {}".format(max_rate))
    print("Fastest rate: {}".format(min_rate))

    spread = max_rate - min_rate
    print("Rate spread: {}\n".format(spread))

    # setting values for partitioning
    upper_value = max_rate
    partitioned_sites_count = 0  # total partitioned sites count
    cutoff = num_chars * 0.1  # cutoff value for creating last partition
    output_phy = "PHYLIP  style"
    output_mrb = "MrBayes style\nbegin mrbayes;"
    output_info = make_output_description(
        divnum, infile_basename, cutoff, max_rate, min_rate, num_chars, spread)
    bin_count = 0
    for partition in range(1, 100):

        Ltest = num_chars - partitioned_sites_count

        if Ltest <= cutoff:  # for last partition to include all the rest
            lower_value = min_rate

            sites = []
            i = 1
            for rate in input_data:
                if upper_value > rate >= lower_value:
                    sites.append(i)
                i += 1

        elif partition == 1:  # for first partition
            lower_value = upper_value - ((upper_value - min_rate) / divnum)
            sites = generate_sites(input_data, lower_value, upper_value)

        else:  # for all other partitions than the last
            lower_value = upper_value - ((upper_value - min_rate) / (divnum + partition * 0.3))
            sites = generate_sites(input_data, lower_value, upper_value)

        bin_count += 1

        partitioned_sites_count += len(sites)  # for total site count

        # info for output in file and screen

        output_info += "\nPartition_{}({} sites):	Rate-span: {}-{}\n".format(
            partition, len(sites), round(upper_value, 6), round(lower_value, 6))  # , BinL

        print("Partition_{} ({} sites): ".format(partition, len(sites)))
        print("Rate-span: {0:.5f}-{0:.5f}\n".format(upper_value, lower_value))

        if sites:
            # setting the output for phylip partitions
            charset = clean_string(str(sites))

            output = "DNA, Partition_{} = {}".format(partition, charset)
            output = clean_string(output)

            output_phy += "\n" + output

            # setting the output format for charsets as MrBayes partitions
            charset = clean_string(str(sites), additional_char=",")

            output = "\nCharset Partition_{} = {};".format(partition, charset)

            output_mrb += clean_string(output, additional_char=",")

        upper_value = lower_value  # resetting the upper range value to the current lower value (for next bin)

        # breaking loop on last partition

        if Ltest <= cutoff:
            break

    # more info
    if partitioned_sites_count != num_chars:
        msg = "Total sites partitioned is not identical to imported " \
              "sites!:{} vs {}".format(partitioned_sites_count, num_chars)
        print(msg)
        output_info += msg

    # fixing partition finishing for output
    # mrb partitioning
    listB = ""
    for partition in range(1, bin_count + 1):
        bapp = "Partition_{}, ".format(partition)
        listB += bapp
    listB = re.sub(", $", "", listB)

    out_finish = "\npartition Partitions = {}: {};".format(bin_count, listB)
    output_mrb += out_finish
    output_mrb += "\nset partition = Partitions;"
    # collecting outputs
    output_finished = [output_info, output_mrb, output_phy]
    output_finished = '\n\n\n'.join(output_finished)
    return output_finished


def generate_sites(input_data, lower_value, upper_value):
    """Generate a list of sites whose evolutionary rate is between upper and lower values"""
    sites = []
    site_index = 1
    for rate in input_data:
        if upper_value >= rate > lower_value:
            sites.append(site_index)
        site_index += 1
    return sites


def make_output_description(divnum, infile_basename, cutoff, max_rate,
                            min_rate, num_chars, spread):
    """Generate text for output for section that describes method"""
    output = "Partition output from ratepartitions.py\n--Written by Tobias Malm " \
         "(20121130)\n\nFor rate file: {} with {} sites!\n".format(
        infile_basename, num_chars)
    output += "\nManually set dividing factor: {}\n".format(divnum)
    output += "Partitions calculated according to:\n\t1st partition: highest " \
          "rate - ((highest rate - minimum_rate)/({0})),\n\tthe remaining as " \
          "lower boundary rate= upper boundary rate -((upper boundary rate - " \
          "minimum rate)/({0}+Partitionnumber*0.3)).\n".format(divnum)
    output += "\tLast partition created when less than 10% of total characters " \
          "are left (={} characters).\n".format(cutoff)
    output += "\nRate spread of entire data set (Highest (slowest, 1=invariant) to " \
          "lowest (fastest) ): Highest: {}, lowest: {}, spread: {}".format(
        max_rate, min_rate, spread)
    return output


def clean_string(text, additional_char=None):
    """Remove unwanted characters from our text"""
    replace_list = ['[', ']', '  ', '\'']

    if additional_char:
        replace_list.append(additional_char)

    for character in replace_list:
        text = text.replace(character, '')
    return text


def read_input_file(infile):
    with open(infile, "r") as handle:
        lines = handle.readlines()

        # reading rates lines into List
        rate_values = []

        for line in lines:
            clean_line = line.strip()
            if clean_line:
                rate_values.append(float(clean_line))

        return rate_values


def write_output_file(output_data, infile, divnum):
    # opening outfile and writing to it
    outfile = "{}_{}.txt".format(infile, divnum)
    print("Output file with rate partition summary and MrBayes and PHYLIP "
          "partition schemes has been created: {}".format(outfile))

    with open(outfile, "w") as handle:
        handle.write(output_data)


def verify_divnum(divnum):
    """Check that dinum is greater or equal than value 1.1"""
    error_msg = "You need to enter factor for division as positive numerical " \
                "value (greater or equal than 1.1)"
    if divnum < 1.1:
        print(error_msg)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "rate_file_txt",
        help="A file with rates separated by hard return (TIGER rate output "
             "-rl command)",
    )
    parser.add_argument(
        "divfactor",
        type=float,
        help="Factor for division as positive numerical value (greater or "
             "equal than 1.1)",
    )

    args = parser.parse_args()
    infile = args.rate_file_txt
    divnum = args.divfactor

    verify_divnum(divnum)

    output_data = run(infile, divnum)
    write_output_file(output_data, infile, divnum)


if __name__ == "__main__":
    main()
