# midas2smol

A code for unpacking MIDAS files directly to SMOL event list files (a minimal-filesize data format used in my [sort codes](https://github.com/e-j-w/GreasySortCodes)). Various bugs are known and the code is an active warzone, user beware.

Based on [grif-replay](https://github.com/GRIFFINCollaboration/grif-replay) by Chris Pearson and Adam Garnsworthy, with my own amateur hour modifications.

## Usage

Compile the code using `make`.


```
midas2smol midas_file output_SMOL_file config_json
```
- `midas_file` should be the path to the first subrun (eg. `run29623_000.mid`).
- `output_SMOL_file` is the event list file to save to disk, which should have the `.smol` file extension.
- `config_json` is the configuration file (in the same JSON format used by grif-replay) which contains pileup correction parameters etc. If no configuration file is specified, the default filename `last.json` is assumed. If the specified configuration file doesn't exist, it will be created.

The data format for output SMOL event lists is specfied in [smol-format.h](smol-format.h).

## Differences from grif-replay

ie. why this is in a separate repo:

- Input and output files are specified using the command line.
- Output is an event list file rather than histograms. The intention is to process these event lists using other [sort codes](https://github.com/e-j-w/GreasySortCodes).
- Addback removed (to be handled by downstream sortcodes).
- grif-angle code removed (to be handled by downstream sortcodes).
- All detector types other than GRIFFIN and BGO were removed (not supported in event list files for now).

## Known issues
- Number of sorted events is not fully consistent each time the code is run.
- Timing is seemingly not consistent with GRSISort, though seems ok for the most part.