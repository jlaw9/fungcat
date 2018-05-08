# fungcat

## Setup on the bioinformatics filesystem
The GO Annotations (GOA) files are located at `/data/inputs/goa`
The STRING networks are located at `/data/inputs/string`
To use those files, add a symbolic link from your inputs folder to those locations:
```
cd inputs
ln -s /data/inputs/goa goa
ln -s /data/inputs/string string
```

The Biorithm/GAIN package can be downloaded from https://bioinformatics.cs.vt.edu/~murali/software/biorithm/index.html
