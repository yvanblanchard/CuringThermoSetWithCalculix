"""
SimpleCure RTM -- Python port of the VBA macro in SimpleCure_RTM.xlsm

Single composite laminate cured in a closed, heated resin-transfer-moulding
tool where BOTH faces are held at the programmed temperature:
  * bottom face  : prescribed temperature = thermal profile
  * top face     : prescribed temperature = thermal profile
  * no surface convection

Note: this workbook uses a shifted Inputs layout -- the fully-cured Tg
threshold is at B19 (row 19) rather than B23 as in the other tools.

Usage:
    python simplecure_rtm.py [workbook.xlsm] [output_dir]
"""

import os
import sys

from simplecure_common import (Inputs, solve_single_domain, write_outputs,
                               print_results)


def run(inp):
    return solve_single_domain(
        inp,
        dirichlet_top=True,         # both faces prescribed = thermal profile
        convection_top=False,       # closed mould, no surface convection
        top_sink_mode=None,
        threshold_row=19,           # fully-cured Tg threshold at B19 (shifted layout)
        series_name_profile="Mould Temperature",
    )


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    workbook = argv[1] if len(argv) > 1 else os.path.join(here, "SimpleCure_RTM.xlsm")
    output_dir = argv[2] if len(argv) > 2 else os.path.join(here, "output_rtm")

    inp = Inputs.from_workbook(workbook)
    out = run(inp)
    print_results("RTM", out["results"])
    write_outputs(out, output_dir)
    print(f"\nField CSVs and chart written to: {output_dir}")


if __name__ == "__main__":
    main(sys.argv)
