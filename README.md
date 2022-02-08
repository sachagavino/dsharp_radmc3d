# dsharp to radmc3d

This routine generates dustkappa*.inp and dustscatmat*.inp files using dsharp_opac package (Birnstiel et al. 2018).


### Parameters:
```ice```: ice fraction should be 0.0 or 0.2.

```amin```: min dust size value in microns.

```amax```: max dust size value in microns.

```nlam```: number of wavelengths. Typically between 200 and 250.

```filename```: choose a filename. Please do not include an underscore (_). Example, if filename = 'myfile', the output filename will be: dustkappa_dsharp-myfile.inp 

```extrapol```: True or False. Extrapolate for large grains to reduce computation time. Default is True.

```scatmat```: True or False. If True, then the matrix elements will be written in the dustkapscatmat_*.inp file. Default value is False. WARNING: don't use True for the moment!

```fig```: yes, no, only. If yes, then a figure is gerenated after the opacities are computed. If only, then the opacity is not computed and only a figure is displayed.

### Example 1:

```
$python dsharp_radmc3d.py -ice 0.2 -amin 5e-3 -amax 1e3 -nlam 250 -filename myfile -scatmat False -fig yes
``` 
--> dust opacities are computed with ice fraction of 0.2, a maximum grain size of 1 mm, 250 wavelengths and a figure of opacity and scattering cross-sections is displayed at the end.

### Example 2:

```
$python dsharp_radmc3d.py -fig only
``` 
--> Opacity tables already exist and the user wants to display the figure only.
