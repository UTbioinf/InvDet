Inversion Detector

version 0.2

In this version, the inverted repeats information is used.

# dependencies

* cmake >= 3.6
* BLASR>=5.2
* MUMmer
* Anything else can be installed via command `pip install -r requirements.txt`

# Install

* go to the root directory of this repository and type `pip install .`
* It's recommended to use virtualenv, e.g., `virtualenv --system-site-packages invdet_env` 
* for private usage, also install aux tools by `./install_aux.sh`, this will install the aux in `build_aux/local/bin/` directory

# Usage:

Type `invdet -h` to see the usage

# TODO

* Extract inverted repeats and view them as validated segments
