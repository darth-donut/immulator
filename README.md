# Immulator

Immunoglobulin simulator.

# Dependencies

1. CMake
2. C/C++ compiler with C++14 support

# Installation

Clone this repository and generate the Makefile using CMake

```bash
$ git clone https://github.com/jfong361/immulator_mkII.git
$ cd immulator_mkii
$ mkdir build && cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```

# Usage

Basic usage, generate 1,000 sequences:

```bash
$ immulator -n 1000 -V human_imgt_ighv -D human_imgt_ighd -J human_imgt_ighj
```

Generate 10,000 sequences and use provided configuration file to simulate the provided abundances of V-(D)-J germline abundances.

```bash
$ immulator -n 10000 -g germ.cfg -V human_imgt_ighv -D human_imgt_ighd -J human_imgt_ighj
```
where germ.cfg should have be a comma separated format as follows:

```bash
$ more germ.cfg
IGHV3-11,70
IGHV2-5,20
IGHV1-18,10
```

which requests that `IGHV3-11` germline gene be simulated in 70% of the sequences (and so on).

## More help

more information about the program can be found using `immulator -h` or `immulator --help`

# Important notes

This repository is very much still in development. Do not expect it to work right out of the box until the first
stable release is announced.
