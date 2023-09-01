# GRFolres

[![status](https://joss.theoj.org/papers/af52e7f1b7637bfa68818fde7c1a34de/status.svg)](https://joss.theoj.org/papers/af52e7f1b7637bfa68818fde7c1a34de)
[![DOI](https://zenodo.org/badge/118786602.svg)](https://zenodo.org/badge/latestdoi/118786602)

GRFolres is an open-source code for performing simulations in modified 
theories of gravity, based on the publicly available 
3+1D numerical relativity code GRChombo.
It is developed and maintained by a collaboration of numerical relativists with a
wide range of research interests, from early universe cosmology to astrophysics
and mathematical general relativity.

GRFolres is written entirely in C++14, using hybrid MPI/OpenMP 
parallelism to achieve good performance on the latest architectures.
It inherits all of the capabilities of the main [GRChombo](https://github.com/GRChombo/GRChombo) 
code, which makes use of the Chombo library for adaptive mesh refinement.

## Getting started
Detailed installation instructions and usage examples are available in
our [wiki](https://github.com/GRChombo/ModifiedGravity/wiki), with the home page giving guidance on where to start.

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the [wiki](https://github.com/GRChombo/ModifiedGravity/wiki)
for our coding style and testing policy before filing a pull request.

## License
GRFolres is licensed under the BSD 3-Clause License. Please see LICENSE for details.

## Citation
A paper for this code is in preparation. In the meantime, please cite the GRChombo JOSS publication using the following bibtex reference:

```
@article{Andrade2021,
  doi = {10.21105/joss.03703},
  url = {https://doi.org/10.21105/joss.03703},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3703},
  author = {Tomas Andrade and Llibert Areste Salo and Josu C. Aurrekoetxea and Jamie Bamber and Katy Clough and Robin Croft and Eloy de Jong and Amelia Drew and Alejandro Duran and Pedro G. Ferreira and Pau Figueras and Hal Finkel and Tiago Fran\c{c}a and Bo-Xuan Ge and Chenxia Gu and Thomas Helfer and Juha Jäykkä and Cristian Joana and Markus Kunesch and Kacper Kornet and Eugene A. Lim and Francesco Muia and Zainab Nazari and Miren Radia and Justin Ripley and Paul Shellard and Ulrich Sperhake and Dina Traykova and Saran Tunyasuvunakool and Zipeng Wang and James Y. Widdicombe and Kaze Wong},
  title = {GRChombo: An adaptable numerical relativity code for fundamental physics},
  journal = {Journal of Open Source Software}
}
```
