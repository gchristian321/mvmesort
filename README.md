Just some basic notes for now. Make this more user friendly later.

To compile, type `make` in the present directory.

To run:
```
./mvmesort <input file> <output file> <channel map> [--no-save-raw] [--no-channel-map-warn]

OPTIONAL ARGUMENTS
  --no--save-raw -->> disable saving of "raw" (detector-level) data to output file (default IS to save)
  --no-channel-map-warn  -->> turn off warnings about problems with the channel map file (default IS to warn)
```

The program will likely only run from the current directory.