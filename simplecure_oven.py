"""
SimpleCure Oven -- Python port of the VBA macro in SimpleCure_Oven.xlsm

Single composite laminate cured on a mould inside a convection oven:
  * bottom face  : prescribed temperature = thermal profile (heated mould)
  * top face     : convection to the oven air (sink = thermal profile), htop

Differs from the Infusion case only in the convection sink temperature (the
oven air follows the programmed thermal profile rather than a fixed ambient),
and this workbook uses cure kinetics model 2 (kin_flag = 2 in the Inputs sheet,
read automatically).

Usage:
    python simplecure_oven.py [workbook.xlsm] [output_dir]
"""

import os
import sys

from simplecure_common import (Inputs, solve_single_domain, write_outputs,
                               print_results)


def run(inp):
    return solve_single_domain(
        inp,
        dirichlet_top=False,        # top face is convective, not prescribed
        convection_top=True,        # convection at the top surface
        top_sink_mode="profile",    # sink temperature = oven air = thermal profile
        threshold_row=23,           # fully-cured Tg threshold at B23
        series_name_profile="Oven Temperature",
    )


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    workbook = argv[1] if len(argv) > 1 else os.path.join(here, "SimpleCure_Oven.xlsm")
    output_dir = argv[2] if len(argv) > 2 else os.path.join(here, "output_oven")

    inp = Inputs.from_workbook(workbook)
    out = run(inp)
    print_results("Oven", out["results"])
    write_outputs(out, output_dir)
    print(f"\nField CSVs and chart written to: {output_dir}")


if __name__ == "__main__":
    main(sys.argv)
