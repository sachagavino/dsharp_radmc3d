# dsharp to radmc3d

This routine generates dustkappa*.inp and dustscatmat*.inp files using dsharp_opac package (Birnstiel et al. 2018).

### Example:

```
$python dsharp_radmc3d.py -ice 0.2 -amin 5e-3 -amax 1e3 -nlam 250 -filename myfile -scatmat False -fig True
``` 

### Parameters:
```ice```: ice fraction should be 0.0 or 0.2.

```amin```: min size value in microns.

```amax```: max size value in microns.

```nlam```: number of wavelengths. Typically between 200 and 250.

```filename```: choose a filename. Please do not include an underscore (_). Example, if filename = 'myfile', the output filename will be: dustkappa_dsharp-myfile.inp 

```scatmat```: True or False. If True, then the matrix elements will be written in the dustkapscatmat_*.inp file. Default value is False.

```fig```: True or False. If True, then a figure is gerenated after the opacities are computed.
