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
    t = 0  # total partitioned sites count
    cutoff = num_chars * 0.10  # cutoff value for creating last partition
    output_phy = ["PHYLIP  style"]
    output_mrb = ["MrBayes style\nbegin mrbayes;"]
    oi = [
        "Partition output from ratepartitions.py\n--Written by Tobias Malm (20121130)\n\nFor rate file: ",
        infile_basename, " with ", str(num_chars), " sites!"]
    oi = ''.join(oi)
    output_info = [oi, '']
    oi = ["Manually set dividing factor: ", str(divnum)]
    oi = ''.join(oi)
    output_info.append(oi)
    oi = [
        "Partitions calculated according to:\n\t1st partition: highest rate - ((highest rate - minimum_rate)/(",
        str(divnum),
        ")),\n\tthe remaining as lower boundary rate= upper boundary rate -((upper boundary rate - minimum rate)/(",
        str(divnum), "+Partitionnumber*0.3)).\n"]
    oi.append(
        "\tLast partition created when less than 10% of total characters are left (=")
    oi.append(str(cutoff))
    oi.append(" characters).")
    oi = ''.join(oi)
    output_info.append(oi)
    output_info.append('')
    oi = [
        "Rate spread of entire data set (Highest (slowest, 1=invariant) to lowest (fastest) ): Highest: ",
        str(max_rate), ", lowest: ", str(min_rate), ", spread: ", str(spread)]
    oi = ''.join(oi)
    output_info.append(oi)
    nBins = 0
    for b in range(1, 100):

        Ltest = num_chars - t

        if (Ltest <= cutoff):  # for last partition to include all the rest

            lowerVal = min_rate

            BinL = []

            i = 1

            for n in input_data:

                if (upper_value > n >= lowerVal):
                    BinL.append(i)

                i += 1

            nBins += 1

        elif (b == 1):  # for first partition

            lowerVal = upper_value - ((upper_value - min_rate) / (divnum))

            BinL = []

            i = 1

            for n in input_data:

                if (upper_value >= n > lowerVal):
                    BinL.append(i)

                i += 1

            nBins += 1

        else:  # for all other partitions than the last

            lowerVal = upper_value - ((upper_value - min_rate) / (divnum + b * 0.3))

            BinL = []

            i = 1

            for n in input_data:

                if (upper_value >= n > lowerVal):
                    BinL.append(i)

                i += 1

            nBins += 1

        t += len(BinL)  # for total site count

        # info for output in file and screen

        oi = ["Partition_", str(b), "(", str(len(BinL)),
              " sites):	Rate-span: ", str(round(upper_value, 6)), "-",
              str(round(lowerVal, 6))]  # , BinL

        oi = ''.join(oi)

        output_info.append(oi)

        output_info.append('')

        pout = ["Partition_", str(b), " (", str(len(BinL)), " sites): "]

        pout = ''.join(pout)

        print(pout)

        pout = ["Rate-span: ", str("{0:.5f}".format(upper_value)), "-",
                str("{0:.5f}".format(lowerVal)), "\n"]

        pout = ''.join(pout)

        print(pout)

        if (len(BinL) > 0):

            # setting the output for phylip partitions

            charset = str(BinL)

            replace_list = ['[', ']', '  ', '\'']

            for rl in replace_list:
                charset = charset.replace(rl, '')

            output = ["DNA, Partition_", str(b), " = ", charset]

            output = ''.join(output)

            for rl in replace_list:
                output = output.replace(rl, '')

            output_phy.append(output)

            # setting the output format for charsets as MrBayes partitions

            charset = str(BinL)

            replace_list = [',', '[', ']', '  ', '\'']

            for rl in replace_list:
                charset = charset.replace(rl, '')

            output = ["Charset Partition_", str(b), " = ", charset, ";"]

            output = ''.join(output)

            for rl in replace_list:
                output = output.replace(rl, '')

            output_mrb.append(output)

        upper_value = lowerVal  # resetting the upper range value to the current lower value (for next bin)

        # breaking loop on last partition

        if (Ltest <= cutoff):
            break

    # more info
    if (t != num_chars):
        print("Total sites paritioned is not identical to imported sites!:", t, " vs ", num_chars)

        oi = ["Total sites paritioned is not identical to imported sites!:",
              str(t), " vs ", str(num_chars)]

        oi = ''.join(oi)

        output_info.append(oi)

    # fixing partition finishing for output
    # mrb partitioning
    listB = []
    for b in range(1, nBins + 1):
        bapp = ["Partition_", str(b)]

        bapp = ''.join(bapp)

        listB.append(bapp)
    listB = ', '.join(listB)
    out_finish = ["partition Partitions = ", str(nBins), ": ", listB, ";"]
    out_finish = ''.join(out_finish)
    output_mrb.append(out_finish)
    output_mrb.append("set partition = Partitions;")
    # collecting outputs
    output_fin1 = '\n'.join(output_info)
    output_fin2 = '\n'.join(output_mrb)
    output_fin3 = '\n'.join(output_phy)
    output_finished = [output_fin1, output_fin2, output_fin3]
    output_finished = '\n\n\n'.join(output_finished)
    return output_finished


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
