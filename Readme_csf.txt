
Cosmological Scaling Function Model (CSF)

First, clone it
$ git clone https://github.com/rodriguezmeza/class_csf.git

Go to CLASS_CSF folder
$ make clean; make all

Then test it:
$ cd tests

Flat
$ time mpirun -n 4 cobaya-run --resume bao_dr2_no-nu_csf_flat.yaml 
$ python getdist_plots_no-nu_csf_flat.py

No-flat
$ time mpirun -n 4 cobaya-run --resume bao_dr2_no-nu_csf_no-flat.yaml 
$ python getdist_plots_no-nu_csf_no-flat.py

