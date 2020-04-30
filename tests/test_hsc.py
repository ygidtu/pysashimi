
import os
from subprocess import check_call, CalledProcessError
import re

from collections import OrderedDict


def __extract_information__(line: str) -> OrderedDict:
    u"""
    extract information from gff of gtf file
    :param line: string after column 8 of gff|gtf file, eg: ID=xxx;Name=xxx or gene_id "xxx"; gene_name "xxx"
    :return: dict
    """
    data = OrderedDict({})
    for i in line.split(";"):
        i = i.strip()
        if not i:
            continue

        tmp = re.split(r"[= ]", i)

        tmp_info = tmp[1].strip() if ":" not in tmp[1] else tmp[1].split(":")[1].strip()
        data[tmp[0]] = tmp_info.replace("\"", "")

    return data


def main(gff, required, gtf):
    data = {}
    with open(gtf) as r:
        for line in r:
            if line.startswith("#"):
                continue

            lines = line.split()

            if lines[2] != "transcript":
                continue

            info = __extract_information__(" ".join(lines[8:]))

            data[info["transcript_id"]] = f"{lines[0]}:{lines[3]}-{lines[4]}:{lines[6]}"

    with open(required) as r:
        required = r.readlines()

    # required = required[5: 15]
    required = [
        # "ENSMUST00000001339",
        "ENSMUST00000001724",
        # "ENSMUST00000023918"
    ]

    for i in set(required):
        i = i.strip()
        try:
            check_call(f"python /mnt/raid61/Personal_data/zhangyiming/code/pysashimi/main.py normal -b /mnt/raid61/PacBio/stats/stringtie/bam_list.txt -g {gff} -o /mnt/raid61/PacBio/stats/stringtie/pdf/8W/{i}.pdf -e {data[i]} -p 20 --color-factor 3 -t 5", shell=True)
        except (KeyError, CalledProcessError):
            continue


if __name__ == '__main__':
    #
    main(
        "/mnt/raid61/PacBio/stats/stringtie/8W_convert.gtf",
        "/mnt/raid61/PacBio/stats/required",
        "/mnt/raid61/PacBio/genome/Mus_musculus/Mus_musculus.GRCm38.93.gtf"
    )