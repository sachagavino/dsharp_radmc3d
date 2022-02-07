# dsharp to radmc3d

This routine generates dustkappa*.inp and dustscatmat*.inp files from dsharp_opac package (Birnstiel et al. 2018).

### Example:

```
$python dsharp_radmc3d.py -ice 0.2 -amin 5e-3 -amax 1e3 -nlam 250 -filename myfile -scatmat True -fig True
``` 
