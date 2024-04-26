## Table of Contents
1. [Installation](#installation)
    1. [Conda](#conda)
    2. [AlleleFinder](#allelefinder)
2. [Tests](#tests)

## Installation <a name="installation"></a>

AlleleFinder is available as a conda package. Therefore, conda must be installed on your system.

### Conda <a name="conda"></a>

Skip this step if you have already installed conda.

\```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
conda update -q conda
\```

### AlleleFinder <a name="allelefinder"></a>

You can now install the AlleleFinder package:

I prefer to create a new conda environment and install the AlleleFinder in a single step:

`conda create -n allelefinder -c olc-bioinformatics allelefinder`

If you wish to create an environment separately:

`conda create -n allelefinder python=3.9`

You can now install AlleleFinder into this environment:

```bash
conda activate allelefinder
conda install -c olcbioinformatics allelefinder
```

If you don't have the ability to install conda for whatever reason, you can install the software with pip.

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

`pip install allelefinder`

You'll also need to install the dependencies. You can do this using the requirements.txt file:

`pip install -r requirements.txt`

To check the installed version of AlleleFinder, use the following command:

`pip show allelefinder | grep Version`

## Tests <a name="tests"></a>

If you encounter issues with the AlleleFinder package, tests are available to ensure that the installation was successful and the basic functionality is present.

You will need to clone this repository and run the tests with pytest. Note: for whatever reason, you need to include the `-s` flag, or you will get an `INTERNALERROR> OSError: [Errno 9] Bad file descriptor` error 

## Running Tests

To ensure that AlleleFinder is working correctly, you can run a set of tests. Here's how to do it:

1. Clone the AlleleFinder repository from GitHub. This will create a copy of the repository on your local machine.

    ```bash
    git clone https://github.com/OLC-Bioinformatics/AlleleFinder.git
    ```

2. Navigate into the cloned repository. This is where the test files are located.

    ```bash
    cd AlleleFinder
    ```

3. Run the tests using pytest. The `-m` option tells Python to run the library module as a script, invoking the module's `__main__` function.

    The `--cov=allele_tools/` option is used to measure code coverage of the `allele_tools` module.

    The `--cov-config=.coveragerc` option specifies the configuration file for coverage.

    The `--cov-report term-missing` option generates a report in the terminal showing lines missed.

    The `-s` option disables all capturing of stdout and stderr by pytest. This means that you can see the print statements from the tests in the console.

    The `-vvv` option increases verbosity. Pytest will give more detailed output, including all test names and the results for each test.

    ```bash
    python -m pytest tests/ --cov=allele_tools/ --cov-config=.coveragerc --cov-report term-missing -s -vvv
    ```