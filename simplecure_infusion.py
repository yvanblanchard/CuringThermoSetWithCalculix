"""
SimpleCure Infusion -- Python port of the VBA macro in SimpleCure_Infusion.xlsm

Single composite laminate cured on a heated mould with resin infusion:
  * bottom face  : prescribed temperature = thermal profile (heated mould)
  * top face     : natural convection to ambient air Too (coefficient htop)

The heavy lifting lives in simplecure_common.solve_single_domain; this script
just supplies the Infusion boundary conditions and I/O.

Usage:
    python simplecure_infusion.py [workbook.xlsm] [output_dir]
"""

import os
import sys

from simplecure_common import (Inputs, solve_single_domain, write_outputs,
                               print_results)


def run(inp):
    return solve_single_domain(
        inp,
        dirichlet_top=False,        # top face is convective, not prescribed
        convection_top=True,        # natural convection at the top surface
        top_sink_mode="ambient",    # sink temperature = Too (Cells(21, 2))
        threshold_row=23,           # fully-cured Tg threshold at B23
        series_name_profile="Mould Temperature",
    )


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    workbook = argv[1] if len(argv) > 1 else os.path.join(here, "SimpleCure_Infusion.xlsm")
    output_dir = argv[2] if len(argv) > 2 else os.path.join(here, "output_infusion")

    inp = Inputs.from_workbook(workbook)
    out = run(inp)
    print_results("Infusion", out["results"])
    write_outputs(out, output_dir)
    print(f"\nField CSVs and chart written to: {output_dir}")


if __name__ == "__main__":
    main(sys.argv)
