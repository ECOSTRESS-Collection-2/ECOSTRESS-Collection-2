.ONESHELL:
SHELL=bash
VERSION := $(shell cat ECOSTRESS/PGEVersion.txt)

# Check whether we have julia available. If so use that julia
ifeq ($(shell conda env list 2> /dev/null | grep ECOSTRESS),ECOSTRESS)
ifeq ($(word 1,$(shell conda list -n ECOSTRESS julia 2> /dev/null | grep julia)),julia)
    JULIA_EXE := $(shell conda run -n ECOSTRESS which julia)
endif
endif
ifeq (${JULIA_EXE},)
    JULIA_EXE := $(shell which julia 2> /dev/null)
endif

default:
	make install

version:
	$(info ECOSTRESS Collection 2 pipeline version ${VERSION})

mamba:
ifeq ($(word 1,$(shell conda run -n base conda list mamba | grep mamba)),mamba)
	@echo "mamba already installed"
else
	-conda deactivate; conda install -y -c conda-forge "mamba>=0.23"
endif

mamba-docker:
ifeq ($(word 1,$(shell conda list mamba | grep mamba)),mamba)
	@echo "mamba already installed"
else
	-conda install -y -c conda-forge mamba
endif

update-env-mamba:
	if ! conda env list | grep '^ECOSTRESS.*ECOSTRESS$$' 2>&1 > /dev/null ; then \
        mamba create -n ECOSTRESS python=3.11 -y; \
    fi; \
    mamba env update -n ECOSTRESS -f ECOSTRESS.yml

# Set up julia for local development
julia-dev-environment:
	${JULIA_EXE} -e 'using Pkg; Pkg.activate("./julia_env"); Pkg.instantiate()'

julia-conda:
	conda install -n ECOSTRESS "julia=1.10" -y
	conda run -n ECOSTRESS curl "https://raw.githubusercontent.com/JuliaLang/MbedTLS.jl/master/src/cacert.pem" --output "${CONDA_PREFIX}/share/julia/cert.pem"

# Basic environment for either development or installation
environment:
	make mamba
	make update-env-mamba
ifeq (${JULIA_EXE},)
	make julia-conda
endif

# Extra packages and preparation for development
dev-environment:
	make environment
	conda run -n ECOSTRESS pip install nose
	make julia-dev-environment

refresh-env:
	make remove
	make environment

environment-docker:
	$(info installing pybind11_cmake inside Docker)
	pip install pybind11_cmake
	$(info installing conda packages inside Docker)
	mamba env update -n base --file ECOSTRESS_docker.yml

clean:
	$(info cleaning build)
	-rm -rvf build
	-rm -rvf dist
	-rm -rvf *.egg-info
	-rm -rvf CMakeFiles
	-rm CMakeCache.txt

uninstall:
	$(info uninstalling ECOSTRESS package)
	-conda run -n ECOSTRESS pip uninstall ECOSTRESS -y

unit-tests:
	$(info running unit tests)
	conda run -n ECOSTRESS nosetests -v -w tests

unit-tests-docker:
	$(info running unit tests inside Docker)
	nosetests -v -w tests

setuptools:
	conda run -n ECOSTRESS python setup.py install

install-julia:
    # This installs into the global base environment
	${JULIA_EXE} -e 'using Pkg; \
	  Pkg.activate(); \
	  Pkg.develop([(;path="./VNP43NRT_jl"), (;path="./STARS_jl")]); \
	  Pkg.add([ \
	    PackageSpec(name="ArchGDAL") \
	    PackageSpec(name="Glob"), \
	    PackageSpec(name="Plots", version="1.40.4"), \
	    PackageSpec(name="Rasters"), \
      ]); \
	  Pkg.activate()'

install-package:
	$(info installing ECOSTRESS package)
	make setuptools
	make clean
	$(info instantiating julia environment)
	make install-julia
	make unit-tests

install-package-docker:
	$(info installing ECOSTRESS package inside Docker)
	python setup.py install
	make clean
	make unit-tests-docker

install:
	make environment
	make clean
	make uninstall
	make install-package

install-docker:
	make clean
	make install-package-docker

# Running this in the base env ensures that we don't fail deletion due to currently
#  executing from the ECOSTRESS environment
remove:
	conda run -n base conda env remove -n ECOSTRESS -y

reinstall-hard:
	make remove
	make install

reinstall-soft:
	make uninstall
	make install-package

docker-build:
	docker build -t pge-eco-level-2-3-4 .

docker-build-mamba:
	docker build --target mamba -t pge-eco-level-2-3-4 .

docker-build-environment:
	docker build --target environment -t pge-eco-level-2-3-4 .

docker-shell:
	docker run -it pge-eco-level-2-3-4 bash

tag-docker-develop-local:
	docker tag pge-eco-level-2-3-4 pge-eco-level-2-3-4:${VERSION}
	docker tag pge-eco-level-2-3-4 cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4
	docker tag pge-eco-level-2-3-4 cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4:${VERSION}

login-docker-develop-local:
	docker login cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

push-docker-develop-local:
	docker push cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4:${VERSION}
	docker push cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

pull-docker-develop-local:
	docker pull cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4
	docker tag cae-artifactory.jpl.nasa.gov:16001/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4 pge-eco-level-2-3-4

docker-develop-local:
	make login-docker-develop-local
	make docker-build
	make tag-docker-develop-local
	make push-docker-develop-local

tag-docker-stage-local:
	docker tag pge-eco-level-2-3-4 cae-artifactory.jpl.nasa.gov:16002/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

login-docker-stage-local:
	docker login cae-artifactory.jpl.nasa.gov:16002/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

push-docker-stage-local:
	docker push cae-artifactory.jpl.nasa.gov:16002/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

pull-docker-stage-local:
	docker pull cae-artifactory.jpl.nasa.gov:16002/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4
	docker tag cae-artifactory.jpl.nasa.gov:16002/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4 pge-eco-level-2-3-4

docker-stage-local:
	make login-docker-stage-local
	make docker-build
	make tag-docker-stage-local
	make push-docker-stage-local

tag-docker-release-local:
	docker tag pge-eco-level-2-3-4 cae-artifactory.jpl.nasa.gov:16003/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

login-docker-release-local:
	docker login cae-artifactory.jpl.nasa.gov:16003/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

push-docker-release-local:
	docker push cae-artifactory.jpl.nasa.gov:16003/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4

pull-docker-release-local:
	docker pull cae-artifactory.jpl.nasa.gov:16003/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4
	docker tag cae-artifactory.jpl.nasa.gov:16003/gov/nasa/jpl/ecostress/sds/pge/pge-eco-level-2-3-4 pge-eco-level-2-3-4

docker-release-local:
	make login-docker-release-local
	make docker-build
	make tag-docker-release-local
	make push-docker-release-local