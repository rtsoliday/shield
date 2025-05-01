# shield

Shield is a C program for performing radiation shielding analyses around a high-energy accelerator. It is a C port of the SHIELD11 program originally developed at SLAC.

## Features

- Defines targets and shielding (primary and secondary) via namelist input
- Computes dose distributions (neutrons and gamma) over a range of angles
- Outputs results in SDDS format for easy plotting and analysis

## Prerequisites

- SDDS library (https://github.com/rtsoliday/SDDS)

## Building

Clone this repository alongside the SDDS repository and ensure the directory structure is:

```
<root>/shield    # this project
<root>/SDDS      # SDDS library sources
```

Then, from the top-level directory:

```bash
make            # builds all dependencies and the shield executable
```

The `shield` binary will be placed in `bin/<OS>-<ARCH>/shield`.

## Usage

```bash
shield <inputFile.in> -dataFile=<material_data.sdds> [options]
```

Options:
- `-describeInput` : Print namelist fields and exit

### Namelist Groups

- `target`           : Define material, dimensions, and attenuation
- `primary_shield`   : Define primary shield material, thickness, position, and angle
- `secondary_shield` : Define additional shielding layers (repeatable)
- `run`              : Beam energy, angle range, number of points, and output file

## Examples

- **Simple run**: `examples/simple1/` demonstrates a single shielding configuration and plotting with `sddsplot`.
- **Scan study**: `examples/scanSecondaryShield1/` shows how to loop over secondary shield thickness.

```bash
# Simple example
cd examples/simple1
env/shield shield.in -dataFile=shieldData.sdds
```

## Authors & Maintenance

- Original program by Michael Borland and Argonne National Laboratory (2016)
- Packaging and maintenance by Robert Soliday
