# How biotic interactions structure species’ responses to perturbations 

[![DOI](https://zenodo.org/badge/737065348.svg)](https://zenodo.org/doi/10.5281/zenodo.10619604)

GitHub repository associated to the [Lajaaiti et al., 2023](change-url) article.

> [!NOTE]
> If you have any questions, please [open an issue](https://github.com/ismael-lajaaiti/species-reactivity/issues) or contact me by e-mail at ismael.lajaaiti@gmail.com.

## Abstract

Predicting how ecological communities will respond to disturbances is a major challenge in community ecology, especially given the variability in species' responses within the same community.
Focusing solely on aggregate responses may obscure extinction risks for certain species due to compensatory effects, emphasizing the need to understand the drivers of the response variability at the species-level. 
Yet, these drivers remain poorly understood. 
Here, we reveal that despite the complexity of the network of biotic interactions, species' responses follow a discernible pattern.
Specifically, we demonstrate that species whose abundance are most reduced by biotic interactions -- which are not always the rarest species -- will be those that exhibit the strongest responses to disturbances. 
This insight enables us to pinpoint sensitive species within communities without requiring precise information about biotic interactions. 
Our novel approach introduces avenues for future research aimed at identifying sensitive species and elucidating their impacts on entire communities.

![pictures/visual-abstract.png](pictures/visual-abstract.png)

## Reproduce figures of the article

To reproduce figures of the article, first you have to clone this repository.

```bash
git clone git@github.com:ismael-lajaaiti/species-reactivity.git
```

Secondly, you can move to the `scripts/` directory and install the required package automatically with:

```bash
cd scripts # From the root of the project.
julia --project=. -e 'using Pkg; Pkg.instantiate()' # Install Julia packages.
```

> [!TIP]
> It is advised to use Julia 1.10, otherwise you may need to resolve conflicts between package versions with:
> ```bash
> julia --project=. -e 'using Pkg; Pkg.update()'
> ```

Finally, you can execute the scripts.

```bash
julia --project=. 01_reactivity-yield-predictability.jl
julia --project=. 02_selection-effect.jl
julia --project=. 03_data.jl
```

Figures are saved in the `scripts/figures/` directory.


## Structure of the code

This project is structured as a Julia package named `SpeciesReactivity`.
The package contains general functions used by the different scripts producing the figures of the article.
For instance, the `response` function allows simulating the response of a community to a specified perturbation.
The code of the package is stored in the `src/` directory.
The code of the package is tested in the `test/` directory.
The scripts producing the figures of the article are stored in the `scripts/` directory.
The figures produced by the scripts are stored in `scripts/figures/` directory (which is created after the first script is run).
