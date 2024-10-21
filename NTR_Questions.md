
Gregory H. Halverson (they/them)
gregory.h.halverson@jpl.nasa.gov
Lead developer
NASA Jet Propulsion Laboratory 329G

Evan W. Davis (he/him)
evan.w.davis@jpl.nasa.gov
Code maintenance 
NASA Jet Propulsion Laboratory 397K

Claire S. Villanueva-Weeks (she/her)
claire.s.villanueva-weeks@jpl.nasa.gov
Code maintenance
NASA Jet Propulsion Laboratory 329G

#### Abstract
This software package includes the code for the ECOSTRESS Collection 2 Jet Propulsion Laboratory Evapotranspiration (JET) suite Product Generating Software (PGE). 

The software was developed as part of a research grant by the NASA Research Opportunities in Space and Earth Sciences (ROSES) program. It was designed for use by the ECOSTRESS mission as a precursor for the Surface Biology and Geology (SBG) mission. However, it may also be useful for general remote sensing and GIS projects in Python. This package can be utilized for remote sensing research in Jupyter notebooks and deployed for operations in data processing pipelines. This software is being released according to the SPD-41 open-science requirements of NASA-funded ROSES projects.

#### This software accomplishes the following:
This software package incldues the glue code and algorithms used for producing the gridded and tiled ECOSTRESS collection 2 data products. This software package produces gridded (HDF5) and tiled (GeoTIFF) outputs for the following data products:
- L2 LSTE: Land Surface Temperature & Emissivity (LSTE)
- L2 STARS: Normalized Difference Vegetation Index (NDVI) & Albedo 
- L3 JET: Evapotranspiration (ET) & Evapotranspiration Auxiliary products (ET_AUX: MET, SM, SEB) 
- L4 JET: Evaporative Stress Index (ESI) & Water Use Efficiency (WUE)

#### What are the unique features of the software?
This is an implemetation of four evapotranpiration algorthims to produce a daily ensemble estimate of evapotranspiration, evaporative stress index, and water use efficiency.

#### What improvements have been made over existing similar software application?
- This is the first production deployment of the JET algorithm suite 
- This software has been optimized for large-scale use for continual processing 

#### What problems are you trying to solve in the software?
This software makes the JET algorithm suite and its related glue code accessible to remote sensing researchers as part of the evapotranspiration modeling capabilities being developed for the ECOSTRESS and SBG missions.

#### Does your work relate to current or future NASA (include reimbursable) work that has value to the conduct of aeronautical and space activities?  If so, please explain:
This software package was developed as part of a research grant by the NASA Research Opportunities in Space and Earth Sciences (ROSES) program. This software was designed for use by the ECOSTRESS mission as a precursor for the SBG mission, but it may be useful generally for remote sensing and GIS projects in python.

#### What advantages does this software have over existing software?
This software can be utilized for remote sensing research in Jupyter notebooks and deployed for operations in data processing pipelines. This is the first software package release of the JET algorithm suite. 

#### Are there any known commercial applications? What are they? What else is currently on the market that is similar?
This software is useful for both remote sensing data analysis and building remote sensing data pipelines.

#### Is anyone interested in the software? Who? Please list organization names and contact information.
NASA ROSES
ECOSTRESS
SBG

#### What are the current hardware and operating system requirements to run the software? (Platform, RAM requirement, special equipment, etc.)
This software is written entirely in python and intended to be distributed using the pip package manager and usable on any system that has access to the following dependencies:
- HDF5
- netCDF4
- GEOS
- GDAL
- TensorFlow

#### How has the software performed in tests? Describe further testing if planned.
This software has been deployed for ECOSTRESS and is currently in production for as a part of the ECOSTRESS data science pipeline.

#### Please identify the customer(s) and sponsors(s) outside of your section that requested and are using your software. 
This package is being released according to the SPD-41 open-science requirements of NASA-funded ROSES projects.

#### Installation
Use the pip package manager to install this package

pip install .
