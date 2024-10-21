# ECOSTRESS Collection 2

This repository contains the code for the ECOSTRESS Collection 2 L2G/L2T-L4G/L4T PGEs.

## Requirements

This system was designed to work in a Linux-like environment and macOS using a conda environment, optionally running within a Docker container.

## Docker

The name of the Docker image is: pge-eco-level-2-3-4
The Artifactory address of the Docker image is: [cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4](cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4)

The 16001 port is for local testing. Use 16002 for integration and testing and 16003 for deployment.

To pull and tag the Docker image:

```bash
$ docker pull cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4
$ docker tag cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4 pge-eco-level-2-3-4
```

The entry-points to execute the PGEs in the Docker image are BASH shell scripts named after each PGE with a .sh extension under a directory called PGE. To run a PGE with Docker:

```bash
$ docker run --rm -it -e PYTHONUNBUFFERED=1 -v /project/sandbox/halverso/ECOSTRESS_15801_013_docker:/working_directory -v /project/sandbox/halverso/ECOSTRESS_15801_013_docker/L1_L2_RAD_LSTE_output:/L1_L2_RAD_LSTE_output pge-eco-level-2-3-4 /bin/bash /pge/L1_L2_RAD_LSTE.sh /working_directory/ECOv002_L1_L2_RAD_LSTE_15801_013_20210419T215859_0700_01_docker.xml
```

To push docker-develop-local on port 16001:

```bash
$ make docker-develop-local
```

To push docker-stage-local on port 16002:

```bash
$ make docker-stage-local
```

To push docker-release-local on port 16003:

```bash
$ make docker-release-local
```


### `conda`

The ECOSTRESS Collection 2 PGEs are designed to run in a Python 3 [`conda`](https://docs.conda.io/en/latest/miniconda.html) environment using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) To use this environment, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Make sure that your shell has been initialized for `conda`.

You should see the base environment name `(base)` when running a shell with conda active.

## Installation

  Use [`git`](https://git-scm.com) to checkout the distribution repository:

```bash
(base) $ git clone git@github.jpl.nasa.gov:halverso/ECOSTRESS-Collection-2.git
```

Navigate to the directory this created:

```bash
(base) $ cd ECOSTRESS-Collection-2
```

Use `make install` to produce the `ECOSTRESS` environment:

```bash
(base) $ make install
```

This should produce a conda environment called `ECOSTRESS` in your [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installation.

## Activation

To use the pipeline, you must activate the `ECOSTRESS` environment:

```bash
(base) $ conda activate ECOSTRESS
```

You should see the environment name `(ECOSTRESS)` in parentheses prepended to the command line prompt.

## Executing PGEs Using Conda Environment

The conda environment includes a command-line entry-point named after each PGE that accepts the filename of the XML run-config as its argument.

```bash
(ECOSTRESS) $ L1_L2_RAD_LSTE run-config.xml
```

The following entry-points are available:
- L1_L2_RAD_LSTE
- L2T_STARS
- L3T_L4T_PTJPL
- L3T_L4T_ALEXI
- L3G_L4G_PTJPL
- L3G_L4G_ALEXI

## Deactivation

When you are done using the pipeline, you can deactivate the `ECOSTRESS` environment:

```bash
(STARS) $ conda deactivate ECOSTRESS
```

You should see the environment name on the command line prompt change to `(base)`.

## Updating

To update your installation of the `ECOSTRESS` environment, first update your copy of the distribution repository:

```bash
(base) $ git pull
```

Then rebuild the `ECOSTRESS` package:

```bash
(base) $ make reinstall-hard
```

## Uninstallation

```bash
(base) $ make remove
```

