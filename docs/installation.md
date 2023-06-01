## Installation

AlleleFinder is available as a conda package, so conda must be installed on your system.

### Conda

Skip this step if you have already installed conda

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
conda update -q conda
```

### AlleleFinder

You can now install the AlleleFinder package:

I prefer to create a new conda environment and install the AlleleFinder in a single step:

`conda create -n allelefinder -c olc-bioinformatics allelefinder=0.1.2=py_0`

If you wish to create an environment separately:

`conda create -n allelefinder python=3.9`

You can now install AlleleFinder into this environment:

```
conda activate allelefinder
conda install -c olcbioinformatics allelefinder=0.1.2=py_0
```

Theoretically, if you don't have the ability to install conda for whatever reason, you can install the software with pip

You'll need a virtual environment with Python 3.9 (start at the appropriate place in these directions depending on your setup):

(if you don't have Python)

`sudo apt install python3`

(if you don't have pip)

`sudo apt install python3-pip`

(if you don't have virtualenv)

`pip install virtualenv`

Create a Python 3.9 virtual environment:

`python3.9 -m venv allelefinder`

Activate the environment:

`source allelefinder/bin/activate`

Install AlleleFinder

`pip install allelefinder==0.1.2`

You'll also need to install the dependencies

```
pip install coloredlogs==15.0.1
pip install geneseekr==0.5.0
```

If you want to run tests, you'll also need to install the testing packages:

```
pip install pytest==7.2.1
pip install pytest-cov==4.0.0
```


## Tests

If you encounter issues with the AlleleFinder package, tests are available to ensure that the installation was successful and the basic functionality is present.

You will need to clone this repository and run the tests with pytest. Note: for whatever reason, you need to include the `-s` flag, or you will get an `INTERNALERROR> OSError: [Errno 9] Bad file descriptor` error 

```
git clone https://github.com/OLC-Bioinformatics/AlleleFinder.git
cd AlleleFinder
python -m pytest tests/ -s
```



