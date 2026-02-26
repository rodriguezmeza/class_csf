
Cosmological Scaling Function Model (CSF)

First, clone it
$ git clone https://github.com/rodriguezmeza/class_csf.git

Go to CLASS_CSF folder
$ make clean; make
$ python setup.py build
$ python setup.py install


Tested in a Mac OS M3 Max with gcc-12 and python 3.12 in an environment. See requirements list by CLASS and install them.

If there is problems generating classy module, then try:

$ python setup.py build_ext --inplace


Then test it:
$ cd tests

Set PATH to class_scf at the beginning of the yaml files.

Flat
$ time mpirun -n 4 cobaya-run --resume bao_dr2_no-nu_csf_flat.yaml 
$ python getdist_plots_no-nu_csf_flat.py

No-flat
$ time mpirun -n 4 cobaya-run --resume bao_dr2_no-nu_csf_no-flat.yaml 
$ python getdist_plots_no-nu_csf_no-flat.py

