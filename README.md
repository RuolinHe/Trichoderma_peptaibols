# Trichoderma_peptaibols
Code for "An Enzymatic Assembly Line Generates a Chemical Arsenal for Precision Interfungal Competition"

Some folders are not included in GitHub due to limitation of file size.

You can get them in [Data and code from: An Enzymatic Assembly Line Generates a Chemical Arsenal for Precision Interfungal Competition](https://doi.org/10.5281/zenodo.18118261)

## Usage
First, you should download the contents in GitHub by:
```
git clone git@github.com:RuolinHe/Trichoderma_peptaibols.git
```

Then, download the data in [Zenodo](https://doi.org/10.5281/zenodo.18118261).

Please decompress antismash.zip, blastp.zip, blastdb.zip, BUSCO.zip and put folders in Trichoderma_peptaibols/data/, and decompress output.zip and put folders in Trichoderma_peptaibols/program/

Finally, your folders should be arranged like this:
```
Trichoderma_peptaibols/
├── .gitignore
├── data/
│   ├── Peptaibol synthetases/
│   ├── blastp/
│   ├── blastdb/
│   ├── BUSCO/
│   ├── antismash/
│   ├── genome/
│   └── other_files
├── program/
│   ├── MLGO/
│   ├── output/
│   ├── reference_tree/
│   ├── species_tree/
│   └── other_files
├── tools_1025/
│   ├── file11.ext
│   └── file12.ext
├── LICENSE.txt
├── README.md
└── script.pptx
```
